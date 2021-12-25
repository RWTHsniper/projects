#include <iostream>
#include <ql/quantlib.hpp>
#include <boost/shared_ptr.hpp>

int main(int, char**) {
    std::cout << "Hello, world!\n";
    // Date Stuff
    QuantLib::Calendar calendar = QuantLib::TARGET();
    QuantLib::Date todaysDate(6, QuantLib::December, 2018);
    QuantLib::Date settlementDate(6, QuantLib::December, 2018);
    QuantLib::Date maturity(6, QuantLib::December, 2019);
    QuantLib::DayCounter dayCounter = QuantLib::Actual365Fixed();
    // Initial Values
    QuantLib::Option::Type OptionType(QuantLib::Option::Call);
    QuantLib::Real stock = 100.00;
    QuantLib::Real strike = 90.00;
    QuantLib::Spread dividendYield = 0.00;
    QuantLib::Rate riskFreeRate = 0.05;
    QuantLib::Volatility volatility = 0.20;

    // Creation of Exercise and Payoff objects:
    boost::shared_ptr<QuantLib::Exercise> europeanExercise(new QuantLib::EuropeanExercise(maturity));
    boost::shared_ptr<QuantLib::StrikedTypePayoff> payoff(new QuantLib::PlainVanillaPayoff(OptionType, strike));
    // Create the European Option (Instrument) object:
    QuantLib::VanillaOption europeanOption(payoff, europeanExercise);
     
    // Variable Handles for the Stochastic (BSM) Process declaration (incuding Flat Term Structures):
    QuantLib::Handle<QuantLib::Quote> spotPrice(boost::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(stock)));
    QuantLib::Handle<QuantLib::YieldTermStructure> flatDividendTS(boost::shared_ptr<QuantLib::YieldTermStructure>(
            new QuantLib::FlatForward(
                settlementDate,
                dividendYield,
                dayCounter)));
    QuantLib::Handle<QuantLib::YieldTermStructure> flatTS(boost::shared_ptr<QuantLib::YieldTermStructure>(
            new QuantLib::FlatForward(
                settlementDate,
                riskFreeRate,
                dayCounter)));
    QuantLib::Handle<QuantLib::BlackVolTermStructure> flatVolTS(boost::shared_ptr<QuantLib::BlackVolTermStructure>(
            new QuantLib::BlackConstantVol(
                settlementDate,
                calendar,
                volatility,
                dayCounter)));

    
    // Creation of the Underlying (1-d) Stochastic Process:
    boost::shared_ptr<QuantLib::BlackScholesMertonProcess> bsmProcess(new QuantLib::BlackScholesMertonProcess(
    spotPrice,
    flatDividendTS,
    flatTS,
    flatVolTS));
    // Creation of the Pricing Engine:
    boost::shared_ptr<QuantLib::PricingEngine> Engine(new QuantLib::AnalyticEuropeanEngine(bsmProcess));
    // Set the Pricing Engine:
    europeanOption.setPricingEngine(Engine);
    
    // Output Results:
    std::cout << "Option Start Date = " << todaysDate << std::endl;
    std::cout << "Option Settlement Date = " << settlementDate << std::endl;
    std::cout << "Option Maturity Date = " << maturity << std::endl;
    std::cout << "Option Type = " << OptionType << std::endl;
    std::cout << "Stock Price = " << stock << std::endl;
    std::cout << "Strike Price = " << strike << std::endl;
    std::cout << "Risk Free Rate = " << riskFreeRate << std::endl;
    std::cout << "Dividend Yield = " << dividendYield << std::endl;
    std::cout << "Volatility = " << volatility << std::endl;
    
    std::cout << "European Option Price = " << europeanOption.NPV() << std::endl;
    
    return 0;
 
    
}
