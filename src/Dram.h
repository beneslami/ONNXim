#ifndef DRAM_H
#define DRAM_H
#include <robin_hood.h>
#include <cstdint>
#include <queue>
#include <utility>
#include <filesystem>

#include "Common.h"
#include "ramulator/Ramulator.hpp"
#include "ramulator2.hh"
#include "memory_system.h" // DRAMsim3

namespace fs = std::filesystem;

class Dram {
 public:
  virtual ~Dram() = default;
  virtual bool running() = 0;
  virtual void cycle() = 0;
  virtual bool is_full(uint32_t cid, MemoryAccess* request) = 0;
  virtual void push(uint32_t cid, MemoryAccess* request) = 0;
  virtual bool is_empty(uint32_t cid) = 0;
  virtual MemoryAccess* top(uint32_t cid) = 0;
  virtual void pop(uint32_t cid) = 0;
  uint32_t get_channel_id(MemoryAccess* request);
  virtual void print_stat() {}
 protected:
  SimulationConfig _config;
  uint32_t _n_ch;
  cycle_type _cycles;
};

// Extended interface — only for backends that support power/bandwidth
class DramWithStats {
public:
    virtual void   collectEpochStats()                               = 0;
    virtual int    getNumChannels()                          const   = 0;
    virtual double getChannelEpochPowerMW(int ch)            const   = 0;
    virtual double getBandwidthGBpsPerChannel(int ch)        const   = 0;
    virtual float  getBandwidthUtilizationPerChannel(int ch) const   = 0;
    virtual double getAggregateBandwidthGBps()               const   = 0;
    virtual float  getAggregateBandwidthUtilization()        const   = 0;
    virtual ~DramWithStats()                                   = default;
};

class SimpleDram : public Dram {
 public:
  SimpleDram(SimulationConfig config);
  virtual bool running() override;
  virtual void cycle() override;
  virtual bool is_full(uint32_t cid, MemoryAccess* request) override;
  virtual void push(uint32_t cid, MemoryAccess* request) override;
  virtual bool is_empty(uint32_t cid) override;
  virtual MemoryAccess* top(uint32_t cid) override;
  virtual void pop(uint32_t cid) override;
 private:
  uint32_t _latency;
  double _bandwidth;

  uint64_t _last_finish_cycle;
  std::vector<std::queue<std::pair<addr_type, MemoryAccess*>>> _waiting_queue;
  std::vector<std::queue<MemoryAccess*>> _response_queue;
};

class DramRamulator : public Dram {
 public:
  DramRamulator(SimulationConfig config);
  virtual bool running() override;
  virtual void cycle() override;
  virtual bool is_full(uint32_t cid, MemoryAccess* request) override;
  virtual void push(uint32_t cid, MemoryAccess* request) override;
  virtual bool is_empty(uint32_t cid) override;
  virtual MemoryAccess* top(uint32_t cid) override;
  virtual void pop(uint32_t cid) override;
  virtual void print_stat() override;
 private:
  std::unique_ptr<ram::Ramulator> _mem;
  robin_hood::unordered_flat_map<uint64_t, MemoryAccess*> _waiting_mem_access;
  std::queue<MemoryAccess*> _responses;

  std::vector<uint64_t> _total_processed_requests;
  std::vector<uint64_t> _processed_requests;
};

class DramRamulator2 : public Dram {
 public:
  DramRamulator2(SimulationConfig config);
  virtual bool running() override;
  virtual void cycle() override;
  virtual bool is_full(uint32_t cid, MemoryAccess* request) override;
  virtual void push(uint32_t cid, MemoryAccess* request) override;
  virtual bool is_empty(uint32_t cid) override;
  virtual MemoryAccess* top(uint32_t cid) override;
  virtual void pop(uint32_t cid) override;
  virtual void print_stat() override;
 private:
  std::vector<std::unique_ptr<NDPSim::Ramulator2>> _mem;
  int _tx_ch_log2;
  int _tx_log2;
  int _req_size;
};

struct ChannelStats {
    double bandwidth_gbps       = 0.0;
    float  bandwidth_util_pct   = 0.0f;
    double power_mw             = 0.0;
};

class DramDRAMsim3 : public Dram, public DramWithStats {
 public:
  DramDRAMsim3(SimulationConfig config, uint64_t freq);
  ~DramDRAMsim3();

  virtual bool running() override;
  virtual void cycle() override;
  virtual bool is_full(uint32_t cid, MemoryAccess* request) override;
  virtual void push(uint32_t cid, MemoryAccess* request) override;
  virtual bool is_empty(uint32_t cid) override;
  virtual MemoryAccess* top(uint32_t cid) override;
  virtual void pop(uint32_t cid) override;
  virtual void print_stat() override;


  void collectEpochStats() override;
  int getNumChannels() const override { return _n_ch; }
  double getChannelEpochPowerMW(int ch) const override { return _mem[ch]->GetEpochPowerMW(); }
  double getBandwidthGBpsPerChannel(int ch) const override { return _stats[ch].bandwidth_gbps; }
  float  getBandwidthUtilizationPerChannel(int ch) const override { return _stats[ch].bandwidth_util_pct; }
  double getAggregateBandwidthGBps() const override;
  float getAggregateBandwidthUtilization() const override;

  double getEpochPowerMW() {
    double total = 0.0;
    for (int ch = 0; ch < _n_ch; ch++) {
        total += _mem[ch]->GetEpochPowerMW();
    }
    return total;
  }
  
 private:
  void read_callback(uint32_t cid, uint64_t addr);
  void write_callback(uint32_t cid, uint64_t addr);

  std::vector<std::unique_ptr<dramsim3::MemorySystem>> _mem;
  std::vector<robin_hood::unordered_map<uint64_t, std::queue<MemoryAccess*>>> _pending;
  std::vector<std::queue<MemoryAccess*>> _response_queue;
  std::vector<uint64_t> _total_processed_requests;
  std::vector<uint64_t> _processed_requests; 
  std::vector<cycle_type> _last_epoch_cycle;
  std::vector<ChannelStats> _stats;
  uint64_t _frequency;
  int _peak_bandwidth_gbps_per_channel;
  int _tx_log2;
  int _tx_ch_log2;
  int _req_size;
};
#endif
