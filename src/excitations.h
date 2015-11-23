#ifndef EXCITATIONS_H
#define EXCITATIONS_H

double PeldorCalculator::det_excitation(double const& freq, exp_signal const& signal_exp) const
{
	double detFreqB1 = 0.25 / signal_exp.detPiHalfLength;
	double detFreqB2 = 0.50 / signal_exp.detPiLength;
	double freqEff(0), freqEffs(0), excitation(0);
	if (detFreqB1 == detFreqB2) {
		freqEff = sqrt(pow(signal_exp.detFreq - freq, 2) + pow(detFreqB1, 2));
		excitation = pow((detFreqB1 / freqEff) * sin(2 * PI * freqEff * signal_exp.detPiHalfLength), 5);
	}
	else {
		freqEff = sqrt(pow(signal_exp.detFreq - freq, 2) + pow(detFreqB1, 2));
		freqEffs = sqrt(pow(signal_exp.detFreq - freq, 2) + pow(detFreqB2, 2));
		excitation = (detFreqB1 / freqEff) * sin(2 * PI * freqEff * signal_exp.detPiHalfLength) * pow((detFreqB2 / freqEffs) * sin(0.5 * 2 * PI * freqEffs * signal_exp.detPiLength), 4);
	}
	return fabs(excitation);
};

double PeldorCalculator::pump_excitation(double const& freq, exp_signal const& signal_exp) const
{
	double pumpFreqB1 = 0.5 / signal_exp.pumpPiLength;
	double freqEff(0), freqEffs(0), excitation(0);
	freqEff = sqrt(pow(signal_exp.pumpFreq - freq, 2) + pow(pumpFreqB1, 2));
	excitation = pow((pumpFreqB1 / freqEff) * sin(0.5 * 2 * PI * freqEff * signal_exp.pumpPiLength), 2);
	return fabs(excitation);
};

#endif