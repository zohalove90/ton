/*
    This file is part of TON Blockchain Library.

    TON Blockchain Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    TON Blockchain Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with TON Blockchain Library.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2017-2020 Telegram Systems LLP
*/

#include "external-message.hpp"
#include "vm/boc.h"
#include "block/block-parse.h"
#include "block/block-auto.h"
#include "block/block-db.h"
#include "fabric.h"
#include "td/actor/actor.h"
#include "td/utils/Random.h"
#include "crypto/openssl/rand.hpp"

namespace ton {

namespace validator {
using td::Ref;

ExtMessageQ::ExtMessageQ(td::BufferSlice data, td::Ref<vm::Cell> root, AccountIdPrefixFull addr_prefix)
    : root_(std::move(root)), addr_prefix_(addr_prefix), data_(std::move(data)) {
  hash_ = block::compute_file_hash(data_);
}

td::Result<Ref<ExtMessageQ>> ExtMessageQ::create_ext_message(td::BufferSlice data) {
  if (data.size() > max_ext_msg_size) {
    return td::Status::Error("external message too large, rejecting");
  }
  vm::BagOfCells boc;
  auto res = boc.deserialize(data.as_slice());
  if (res.is_error()) {
    return res.move_as_error();
  }
  if (boc.get_root_count() != 1) {
    return td::Status::Error("external message is not a valid bag of cells");  // not a valid bag-of-Cells
  }
  auto ext_msg = boc.get_root_cell();
  if (ext_msg->get_level() != 0) {
    return td::Status::Error("external message must have zero level");
  }
  if (ext_msg->get_depth() >= max_ext_msg_depth) {
    return td::Status::Error("external message is too deep");
  }
  vm::CellSlice cs{vm::NoVmOrd{}, ext_msg};
  if (cs.prefetch_ulong(2) != 2) {  // ext_in_msg_info$10
    return td::Status::Error("external message must begin with ext_in_msg_info$10");
  }
  ton::Bits256 hash{ext_msg->get_hash().bits()};
  if (!block::gen::t_Message_Any.validate_ref(128, ext_msg)) {
    return td::Status::Error("external message is not a (Message Any) according to automated checks");
  }
  if (!block::tlb::t_Message.validate_ref(128, ext_msg)) {
    return td::Status::Error("external message is not a (Message Any) according to hand-written checks");
  }
  block::gen::CommonMsgInfo::Record_ext_in_msg_info info;
  if (!tlb::unpack_cell_inexact(ext_msg, info)) {
    return td::Status::Error("cannot unpack external message header");
  }
  auto dest_prefix = block::tlb::t_MsgAddressInt.get_prefix(info.dest);
  if (!dest_prefix.is_valid()) {
    return td::Status::Error("destination of an inbound external message is an invalid blockchain address");
  }
  return Ref<ExtMessageQ>{true, std::move(data), std::move(ext_msg), dest_prefix};
}

void ExtMessageQ::run_message(td::BufferSlice data, td::actor::ActorId<ton::validator::ValidatorManager> manager,
                              td::Promise<td::Unit> promise) {
  auto R = create_ext_message(std::move(data));
  if (R.is_error()) {
    return promise.set_error(R.move_as_error_prefix("failed to parse external message "));
  }
  auto M = R.move_as_ok();
  auto root = M->root_cell();
  block::gen::CommonMsgInfo::Record_ext_in_msg_info info;
  tlb::unpack_cell_inexact(root, info); // checked in create message
  ton::StdSmcAddress addr;
  ton::WorkchainId wc;
  if(!block::tlb::t_MsgAddressInt.extract_std_address(info.dest, wc, addr)) {
    return promise.set_error(td::Status::Error(PSLICE() << "Can't parse destination address"));
  }

  run_fetch_account_state(wc, addr, manager,
      [promise = std::move(promise), msg_root = root, wc = wc](td::Result<std::tuple<td::Ref<vm::CellSlice>,UnixTime,LogicalTime,std::unique_ptr<block::ConfigInfo>>> res) mutable {
        if (res.is_error()) {
          promise.set_error(td::Status::Error(PSLICE() << "Failed to get account state"));
        } else {
          auto tuple = res.move_as_ok();
          block::Account acc;
          auto shard_acc = std::get<0>(tuple);
          auto utime = std::get<1>(tuple);
          auto lt = std::get<2>(tuple);
          auto config = std::move(std::get<3>(tuple));
          if(!acc.unpack(shard_acc, {}, utime, false)) {
            promise.set_error(td::Status::Error(PSLICE() << "Failed to unpack account state"));
          }
          if(run_message_on_account(wc, acc, utime, lt + 1, msg_root, std::move(config))) {
            promise.set_value(td::Unit());
          } else {
            promise.set_error(td::Status::Error(PSLICE() << "External message was not accepted"));
          }
        }
      }
  );
}

bool ExtMessageQ::run_message_on_account(ton::WorkchainId wc,
                                         block::Account acc,
                                         UnixTime utime, LogicalTime lt,
                                         td::Ref<vm::Cell> msg_root,
                                         std::unique_ptr<block::ConfigInfo> config) {

   std::vector<block::StoragePrices> storage_prices_;
   td::BitArray<256> rand_seed_;
   block::ComputePhaseConfig compute_phase_cfg_;
   block::ActionPhaseConfig action_phase_cfg_;
   {
     auto res = config->get_storage_prices();
     if (res.is_error()) {
       LOG(DEBUG) << "Can not unpack storage prices";
       return false;
     }
     storage_prices_ = res.move_as_ok();
   }
   block::StoragePhaseConfig storage_phase_cfg_{&storage_prices_};
   {
    // generate rand seed
    prng::rand_gen().strong_rand_bytes(rand_seed_.data(), 32);
    LOG(DEBUG) << "block random seed set to " << rand_seed_.to_hex();
   }
   {
    // compute compute_phase_cfg / storage_phase_cfg
     auto cell = config->get_config_param(wc == ton::masterchainId ? 20 : 21);
     if (cell.is_null()) {
       LOG(DEBUG) << "cannot fetch current gas prices and limits from masterchain configuration";
       return false;
     }
     if (!compute_phase_cfg_.parse_GasLimitsPrices(std::move(cell), storage_phase_cfg_.freeze_due_limit,
                                                  storage_phase_cfg_.delete_due_limit)) {
       LOG(DEBUG) <<"cannot unpack current gas prices and limits from masterchain configuration";
       return false;
     }
     compute_phase_cfg_.block_rand_seed = rand_seed_;
     compute_phase_cfg_.libraries = std::make_unique<vm::Dictionary>(config->get_libraries_root(), 256);
     compute_phase_cfg_.global_config = config->get_root_cell();
   }
   {
    // compute action_phase_cfg
    block::gen::MsgForwardPrices::Record rec;
    auto cell = config->get_config_param(24);
    if (cell.is_null() || !tlb::unpack_cell(std::move(cell), rec)) {
      LOG(DEBUG) << "cannot fetch masterchain message transfer prices from masterchain configuration";
      return false;
    }
    action_phase_cfg_.fwd_mc =
        block::MsgPrices{rec.lump_price,           rec.bit_price,          rec.cell_price, rec.ihr_price_factor,
                         (unsigned)rec.first_frac, (unsigned)rec.next_frac};
    cell = config->get_config_param(25);
    if (cell.is_null() || !tlb::unpack_cell(std::move(cell), rec)) {
      LOG(DEBUG) << "cannot fetch standard message transfer prices from masterchain configuration";
      return false;
    }
    action_phase_cfg_.fwd_std =
        block::MsgPrices{rec.lump_price,           rec.bit_price,          rec.cell_price, rec.ihr_price_factor,
                         (unsigned)rec.first_frac, (unsigned)rec.next_frac};
    action_phase_cfg_.workchains = &config->get_workchain_list();
    action_phase_cfg_.bounce_msg_body = (config->has_capability(ton::capBounceMsgBody) ? 256 : 0);
  }

  std::unique_ptr<block::Transaction> trans =
      std::make_unique<block::Transaction>(acc, block::Transaction::tr_ord, lt, utime, msg_root);
  bool ihr_delivered = false;  // FIXME
  if (!trans->unpack_input_msg(ihr_delivered, &action_phase_cfg_)) {
      // inbound external message was not accepted
      LOG(DEBUG) << "inbound external message rejected by account"
                 << " before smart-contract execution";
      return false;
  }

  if (!trans->prepare_storage_phase(storage_phase_cfg_, true, true)) {
      LOG(DEBUG) << "cannot create storage phase of a new transaction for smart contract ";
      return false;
  }
  if (!trans->prepare_compute_phase(compute_phase_cfg_)) {
    LOG(DEBUG) << "cannot create compute phase of a new transaction for smart contract ";
    return false;
  }
  if (!trans->compute_phase->accepted) {
    // inbound external message was not accepted
    LOG(DEBUG) << "inbound external message rejected by transaction ";
    return false;
  }
  if (trans->compute_phase->success && !trans->prepare_action_phase(action_phase_cfg_)) {
    LOG(DEBUG) << "cannot create action phase of a new transaction for smart contract ";
    return false;
  }
  if (!trans->serialize()) {
    LOG(DEBUG) << "cannot serialize new transaction for smart contract ";
    return false;
  }
  auto trans_root = trans->commit(acc);
  if (trans_root.is_null()) {
    LOG(DEBUG) << "cannot commit new transaction for smart contract ";
    return false;
  }
  return true;

}

}  // namespace validator
}  // namespace ton
