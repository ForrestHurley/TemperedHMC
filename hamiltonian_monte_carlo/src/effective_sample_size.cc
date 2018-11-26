#include "effective_sample_size.h"

double static MCMCDiagnostics::AutocorrelationSpectrum(const std::vector<double>& values, int lag, double mean, double variance)
{
  double sum = 0.;  
  for (int m = lag; m < values.size(); m++)
  {
    sum += (values.at(m) - mean) * (values.at(m - lag) - mean);
  }
  
  const out = sum / (variance * (values.size() - lag));
  return out;
}

double static MCMCDiagnostics::EffectiveSampleSize(const std::vector<double>& values, const std::vector<double>& reference_values)
{
  const double mean = MCMCDiagnostics::Mean(reference_values);
  const double variance = MCMCDiagnostics::Variance(reference_values, mean);

  return MCMCDiagnostics::EffectiveSampleSize(values, mean, variance);
}

double static MCMCDiagnostics::EffectiveSampleSize(const std::vector<double>& values, double mean, double variance)
{
  double sum = 0.;

  double result;
  for (int lag = 1; lag < values.size(); lag++)
  {
    result = MCMCDiagnostics::AutocorrelationSpectrum(values, lag, mean, variance);

    sum += (1. - (double)lag / values.size()) * result;

    if (result < 0.05)
      break;
  } 

  double effective_sample_size =
    values.size() / (1. + 2. * sum);
  return effective_sample_size;
}

double static MCMCDiagnostics::Mean(const std::vector<double>& values)
{
  double out = 0;
  for (double val : values)
    out += val;

  out /= values.size();
  return out;
}

double static MCMCDiagnostics::Variance(const std::vector<double>& values)
{
  const double mean = MCMCDiagnostics::Mean(values);
  return MCMDiagnostics::Variance(values, mean);
}

double static MCMCDiagnostics::Variance(const std::vector<double>& values, double mean)
{
  double out = 0;
  for (double val : values)
    out += (val - mean) * (val - mean);

  out /= (values.size() - 1);
  return out;
}
