#ifndef MCMC_DIAGNOSTICS_H
#define MCMC_DIAGNOSTICS_H

#include <vector>

class MCMCDiagnostics
{
private:
  double static AutocorrelationSpectrum(const std::vector<double>& values, int lag, double mean, double variance);
public:
  double static EffectiveSampleSize(const std::vector<double>& values, const std::vector<double>& reference_values);
  double static EffectiveSampleSize(const std::vector<double>& values, double mean, double variance);

  double static Mean(const std::vector<double>& values);
  double static Variance(const std::vector<double>& values);
  double static Variance(const std::vector<double>& values, double mean);
};

#endif
