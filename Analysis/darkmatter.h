#ifndef analysis_darkmatter_h
#define analysis_darkmatter_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"

namespace MA5
{
class darkmatter : public AnalyzerBase
{
  INIT_ANALYSIS(darkmatter,"darkmatter")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);
  virtual double reliso04(const EventFormat& event, const RecLeptonFormat& lep);

 private:
    std::vector<const RecLeptonFormat*> selectedels;
    std::vector<const RecLeptonFormat*> selectedmus;
    std::vector<const RecJetFormat*> selectedjets;
    int nevent;
    int ee;
    int mm;

    float Zmass = 91.876;
};
}

#endif
