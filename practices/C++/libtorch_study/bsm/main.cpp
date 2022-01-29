/*
Reference sites
https://pytorch.org/tutorials/advanced/cpp_autograd.html
https://qvious.com/docs/1-1_pricing/

*/

#pragma warning(disable:4251 4275 4244 4267 4522 4273 4018 4305)
#include <torch/torch.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <ctime>
#include <chrono>
#include <random>
#define MAX(x,y) (((x)>(y)) ? (x) : (y))

class Timer {
public:
	Timer() : beg_(clock_::now()) {}
	void reset() { beg_ = clock_::now(); }
	double elapsed() const {
		return std::chrono::duration_cast<second_>
			(clock_::now() - beg_).count();
	}
private:
	typedef std::chrono::high_resolution_clock clock_;
	typedef std::chrono::duration<double, std::ratio<1> > second_;
	std::chrono::time_point<clock_> beg_;
};

enum OptionType {
	Call = 1, Put = -1
};

struct PriceResult {
	double price;
	double delta;
	double gamma = 0;
};

double normcdf(double x, double mu=0, double sigma=1) {
	double v = (x - mu) / sigma;
	return 0.5 * erfc(-v * M_SQRT1_2);
}

PriceResult bsprice(double s, double k,	double r, double q, double t,
   double sigma, OptionType type) {
	double d1 = (log(s / k) + (r - q + 0.5 * sigma * sigma) * t) / (sigma * sqrt(t));
	double d2 = d1 - sigma * sqrt(t);
	double nd1 = normcdf(type * d1);
	double nd2 = normcdf(type * d2);
	double price = type * (s * exp(-q * t) * nd1 - k * exp(-r * t) * nd2);
	double delta = type * exp(-q * t) * nd1;
	return PriceResult({ price, delta });
}

PriceResult mcprice_cpu(double s, double k, double r, double q, double t, 
  double sigma, OptionType type, unsigned int ntimes, unsigned int numOfSimulation) {
	double sumOfPayoff = 0, sumOfDelta = 0;
	double df = exp(-r * t);
	std::mt19937_64 gen;
	std::normal_distribution<double> engine(0.0, 1.0);
	gen.seed(std::random_device{}());
	double dt = t / ntimes;
	double es = exp((r - q - 0.5 * sigma * sigma) * dt);
	double diffution = sigma * sqrt(dt);
	for (unsigned int i = 0; i < numOfSimulation; ++i) {
		double st = s;
		for (unsigned int j = 0; j < ntimes; ++j) {
			double e = engine(gen);
			st = st * es * exp(diffution * e);
		}
		double p = MAX(type * (st - k), 0);
		double d = type * (st - k) > 0 ? type * st / s : 0;
		sumOfPayoff += df * p;
		sumOfDelta += df * d;
	}
	double price = sumOfPayoff / numOfSimulation;
	double delta = sumOfDelta / numOfSimulation;
	return PriceResult({ price, delta });
}

PriceResult mcprice_torch(double s0, double k, double r, double q, double t,
  double sigma, OptionType type, unsigned int ntimes, unsigned int numOfSimulation) {
	torch::manual_seed(rand());

	// Create the device we pass around based on whether CUDA is available.
	torch::Device device(torch::kCPU);
	if (torch::cuda::is_available()) {
		//std::cout << "CUDA is available! Training on GPU." << std::endl;
		device = torch::Device(torch::kCUDA);
	}

	auto options_0 =
		torch::TensorOptions()
		.dtype(torch::kFloat64)
		.device(device)
		.requires_grad(true);

	auto options_1 =
		torch::TensorOptions()
		.dtype(torch::kFloat64)
		.device(device)
		.requires_grad(false);

	int m = numOfSimulation;

	torch::Tensor s_0 = torch::tensor(s0, options_0);

	double dt = t / ntimes;
	//int ntimes = int(t / dt);
	torch::Tensor ts = torch::arange(0.0, t + dt, dt, options_1); // first elem:0, last elem: t

	torch::Tensor e = torch::randn({ m, ntimes }, options_1);
	torch::Tensor sum_e = torch::cumsum(e, 1); // [0-th e[0], e[0]+e[1], ..., e[0]+e[n]+...e[m-1]; 1-th ...;...]
	sum_e = torch::cat({ torch::zeros({ m, 1 }, options_1), sum_e }, 1); // concatenate a col of zeros + sum_e
	torch::Tensor s = s_0 * torch::exp((r - q - 0.5 * sigma * sigma) * ts + sigma * sqrt(dt) * sum_e);

	double disc = exp(-r * t);
	auto sAtMat = s.index({ torch::arange(0, s.size(0), torch::kLong), torch::tensor(s.size(1) - 1) }); // last column. s[:,-1]
	// torch::Tensor discPayoff = disc * torch::relu(type *( sAtMat - k)); // original version
	torch::Tensor discPayoff = disc * at::softplus(type *( sAtMat - k), 100.0, 10.0); // Softplus
	auto price = discPayoff.mean();

	// Automatic differentiation
    auto grad_output = torch::ones_like(price); // jacobian
	auto first_derivative = torch::autograd::grad({price}, {s_0}, /*grad_outputs=*/{}, /*retain_graph=*/c10::nullopt, /*create_graph=*/true)[0];
	std::cout << first_derivative.requires_grad() << std::endl;
    auto second_derivative = torch::autograd::grad({first_derivative}, {s_0}, /*grad_outputs=*/{}, /*retain_graph=*/c10::nullopt, /*create_graph=*/true)[0];
	// price.backward();

	// return PriceResult({ price.item<double>(), s_0.grad().item<double>() }); // original version
	return PriceResult({ price.item<double>(), first_derivative.item<double>(), second_derivative.item<double>() });
}

int main() {

	Timer tmr;
	double s = 100, k = 100, r = 0.02, q = 0.01, t = 0.25, sigma = 0.15;
	// s = 90.0;
	OptionType type = Call;
	unsigned int ntimes = 100, numOfSimulation = 50000;

	std::cout << "Pricing with BS Formula" << std::endl;
	tmr.reset();
	PriceResult res_anal = bsprice(s, k, r, q, t, sigma, type);	
	std::cout << "Anal Price = " << res_anal.price << "\t";
	std::cout << "Delta = " << res_anal.delta << std::endl;
    std::cout << "Gamma = " << res_anal.gamma << std::endl;
	std::cout << std::string(40, '-') << std::endl;
		
	tmr.reset();
	std::cout << "Pricing with STD" << std::endl;
	PriceResult res_std = mcprice_cpu(s, k, r, q, t, sigma, type, ntimes, numOfSimulation);
	double computationTime = tmr.elapsed();
	std::cout << "MC Price = " << res_std.price << "\t";

	std::cout << "Delta = " << res_std.delta << std::endl;
    std::cout << "Gamma = " << res_std.gamma << std::endl;
	std::cout << "Time = " << computationTime << std::endl;
	std::cout << std::string(40, '-') << std::endl;

	
	for (int _i = 0; _i < 10; _i++) {
		tmr.reset();		
		std::cout << "Pricing with Torch  #" << _i+1 << std::endl;
		PriceResult res_torch = mcprice_torch(s, k, r, q, t, sigma, type, ntimes, numOfSimulation);
		computationTime = tmr.elapsed();
		std::cout << "MC Price = " << res_torch.price << "\t";
        // std::cout << "Delta = " << res_anal.delta << "\t" << "Gamma = " << res_torch.gamma << std::endl;
		std::cout << "Delta = " << res_torch.delta << std::endl;
		std::cout << "Gamma = " << res_torch.gamma << std::endl;
		std::cout << "Time = " << computationTime << std::endl;
		std::cout << std::string(40, '-') << std::endl;
	}	

  return 0;
}