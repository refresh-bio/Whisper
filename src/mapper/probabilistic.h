// *******************************************************************************************
// This file is a part of Whisper reads mapper.
// The homepage of the project is http://sun.aei.polsl.pl/REFRESH/whisper
// 
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Adam Gudys
// 
// Version : 1.1
// Date    : 2018-07-10
// License : GNU GPL 3
// *******************************************************************************************


#pragma once

// ************************************************************************************
class IDistribution 
{
public:
	virtual double pdf(double x) const = 0;
	virtual double cdf(double x) const = 0;
	virtual double invcdf(double x) const = 0;
};

// ************************************************************************************
// Normal distribution wrapper
class NormalDistribution : public IDistribution 
{
protected:
	double meanX;
	double devX;
	double meanX2;

public:
	virtual double mean() const { return meanX; }
	virtual double dev() const { return devX; }

	/// Initializes mean and standard deviation
	NormalDistribution(double mean, double dev) : meanX(mean), devX(dev), meanX2(mean*mean) {}

	/// Probability density function
	virtual double pdf(double x) const override { return 0; }

	/// Cumulative distribution function (from Handbook of Mathematical Functions)
	virtual double cdf(double x) const override {
		// get standardized argument
		double sx = (x - meanX) / devX;
		double out = 0.5 * std::erfc(-sx * M_SQRT1_2);
		return out;
	}

	/// Inversed cumulative distribution function (from Handbook of Mathematical Functions)
	virtual double invcdf(double p) const override {
		
		if (p <= 0.0 || p >= 1.0) {
			throw std::invalid_argument("Invalid invcdf() argument");
		}

		double X = p < 0.5
			? -rationalApproximation(std::sqrt(-2.0 * std::log(p))) // F^-1(p) = - G^-1(p)
			: rationalApproximation(std::sqrt(-2.0 * std::log(1 - p))); 	// F^-1(p) = G^-1(1-p)
			
		// de-standardize
		X = X * devX + meanX;
		return X;
	}

	void addSample(double x, int n) {
		double coeff = (double)(n) / (n + 1);

		meanX = meanX * coeff + x / (n + 1);
		meanX2 = meanX2 * coeff + x * x / (n + 1);

		devX = std::sqrt(meanX2 - meanX * meanX);
	}

protected:
	double rationalApproximation(double t) const
	{
		// Abramowitz and Stegun formula 26.2.23.
		// The absolute value of the error should be less than 4.5 e-4.
		const double c[] = { 2.515517, 0.802853, 0.010328 };
		const double d[] = { 1.432788, 0.189269, 0.001308 };
		return t - ((c[2] * t + c[1])*t + c[0]) /
			(((d[2] * t + d[1])*t + d[0])*t + 1.0);
	}
};

// ************************************************************************************
class InsertSizeModel 
{
protected:
	scoring_t scoring;
	size_t samplesCount;
	size_t minSamplesCount;
	size_t maxSamplesCount;
	double maxPenalty;
	double penaltySaturationSigmas;
	double highConfidenceSigmas;
	NormalDistribution initialDistribution;
	NormalDistribution currentDistibution;

public:
	InsertSizeModel(
		double mean, double dev, int minSamplesCount, int maxSamplesCount,
		double highConfidenceSigmas, double penaltySaturationSigmas,
		scoring_t scoring)
		: scoring(scoring),
		minSamplesCount(minSamplesCount),
		maxSamplesCount(maxSamplesCount),
		samplesCount(0), 
		penaltySaturationSigmas(penaltySaturationSigmas),
		highConfidenceSigmas(highConfidenceSigmas),
		initialDistribution(mean, dev), 
		currentDistibution(mean, dev)		
	{}

	virtual ~InsertSizeModel()
	{}

	const scoring_t & getScoring() const { return scoring; }

	size_t getMinSamplesCount() const { return minSamplesCount; }

	size_t getSamplesCount() const { return samplesCount; }

	double getHighConfidenceDev() const { return highConfidenceSigmas * dev(); }

	double getMyersDev() const { return penaltySaturationSigmas * dev(); }

	double getSaturationDev() const { return penaltySaturationSigmas * dev(); }

	double getMinPenalty() const { return calculatePenalty(this->mean(), false); }

	double getMaxPenalty() const { return calculatePenalty(this->mean() - getSaturationDev(), false); }

	void reset() {
		samplesCount = 0;
		currentDistibution = initialDistribution;
	}

	/// Use samples to recalculate internal mean;
	void addSample(double x) {
		currentDistibution.addSample(x, samplesCount);
		++samplesCount;
	}

	/// Calculates penalty on the basis of insert size
	double calculatePenalty(int insertSize, bool saturated = true) const {
		
		if (saturated) {
			insertSize = std::min(int(mean() + getSaturationDev()), insertSize);
			insertSize = std::max(int(mean() - getSaturationDev()), insertSize);
		}

		double p = (insertSize < mean()) ?
			cdf((double)insertSize) :
			cdf(mean() - ((double)insertSize - mean()));

		if (p == 0) {
			p = std::numeric_limits<double>::min();
		}

		double log_p = log_base(p, 4.0);
		double insertSizePenalty = -scoring.match * log_p;

		return insertSizePenalty;
	}

	/// Calculates insert size on the basis of penalty
	int calculateInvPenalty(double penalty) {
		int tlen = (int)invcdf(pow(4.0, -penalty / scoring.match));
		return tlen;
	}

	virtual double mean() const {
		return currentDistibution.mean();
	}

	virtual double dev() const {
		return currentDistibution.dev();
	}

	virtual double pdf(double x) const {
		return currentDistibution.pdf(x);
	}

	virtual double cdf(double x) const {
		return currentDistibution.cdf(x);
	}

	virtual double invcdf(double x) const {
		return currentDistibution.invcdf(x);
	}
};

// EOF
