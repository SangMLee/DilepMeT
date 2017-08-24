#ifndef analysis_DilepMet_h
#define analysis_DilepMet_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"

namespace MA5
{
class DilepMet : public AnalyzerBase
{
  INIT_ANALYSIS(DilepMet,"DilepMet")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

 private:
  double myEvent;
};
}

#endif
