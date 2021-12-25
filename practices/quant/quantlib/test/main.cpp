#include <iostream>
#include <ql/quantlib.hpp>


#include <ql/qldefines.hpp>
#if !defined(BOOST_ALL_NO_LIB) && defined(BOOST_MSVC)
#  include <ql/auto_link.hpp>
#endif
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/vanilla/binomialengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
#include <ql/pricingengines/vanilla/baroneadesiwhaleyengine.hpp>
#include <ql/pricingengines/vanilla/bjerksundstenslandengine.hpp>
#include <ql/pricingengines/vanilla/batesengine.hpp>
#include <ql/pricingengines/vanilla/integralengine.hpp>
#include <ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp>
#include <ql/pricingengines/vanilla/mceuropeanengine.hpp>
#include <ql/pricingengines/vanilla/mcamericanengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanvasicekengine.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/models/shortrate/onefactormodels/vasicek.hpp>


namespace ql = QuantLib;
void testSolver1(){
    std::cout << "Example from the pdf" << std::endl;

    // ql::Real S0 =100.0 , K =105.0;
    // ql::Real rd =0.034 , rf =0.021 , tau =0.5 , vol =0.177;
    ql::Real S0 =100.0 , K =90.0;
    ql::Real rd =0.05 , rf =0.0 , tau =1.0 , vol =0.2; // In this case, rf is used as a dividend yield
    ql::Real domDisc = std :: exp (-rd* tau ), forDisc = std :: exp (-rf* tau );
    ql::Real stdDev = vol * std :: sqrt ( tau );
    boost :: shared_ptr < ql::PlainVanillaPayoff > vanillaPayoffCall (
        new ql::PlainVanillaPayoff ( ql::Option :: Call ,K ));
    boost :: shared_ptr < ql::PlainVanillaPayoff > vanillaPayoffPut (
        new ql::PlainVanillaPayoff ( ql::Option :: Put ,K ));
    boost :: shared_ptr < ql::AssetOrNothingPayoff > aonPayoffCall (
        new ql::AssetOrNothingPayoff ( ql::Option :: Call ,K ));
    ql::BlackScholesCalculator vanillaCallPricer ( vanillaPayoffCall ,S0 , forDisc ,stdDev , domDisc );
    ql::BlackScholesCalculator vanillaPutPricer ( vanillaPayoffPut ,S0 , forDisc ,stdDev , domDisc );
    ql::BlackScholesCalculator aonCallPricer ( aonPayoffCall ,S0 , forDisc , stdDev , domDisc );
    std :: cout << " -------------- Vanilla Call Values -------------" << std :: endl ;
    std :: cout << " Value :" << vanillaCallPricer . value () << std :: endl ;
    std :: cout << " Delta :" << vanillaCallPricer . delta () << std :: endl ;
    std :: cout << " Gamma :" << vanillaCallPricer . gamma () << std :: endl ;
    std :: cout << " Vega :" << vanillaCallPricer . vega (tau ) << std :: endl ;
    std :: cout << " Theta :" << vanillaCallPricer . theta ( tau ) << std :: endl ;
    std :: cout << " Delta Fwd :" << vanillaCallPricer . deltaForward () << std :: endl ;
    std :: cout << " Gamma Fwd :" << vanillaCallPricer . gammaForward () << std :: endl ;
    std :: cout << " -------------- Vanilla Put Values -------------" << std :: endl ;
    std :: cout << " Value :" << vanillaPutPricer . value () << std :: endl ;
    std :: cout << " Delta :" << vanillaPutPricer . delta () << std :: endl ;
    std :: cout << " Gamma :" << vanillaPutPricer . gamma () << std :: endl ;
    std :: cout << " Vega :" << vanillaPutPricer . vega (tau ) << std :: endl ;
    std :: cout << " Theta :" << vanillaPutPricer . theta ( tau ) << std :: endl ;
    std :: cout << " Delta Fwd :" << vanillaPutPricer . deltaForward () << std :: endl ;
    std :: cout << " Gamma Fwd :" << vanillaPutPricer . gammaForward () << std :: endl ;
    std :: cout << " -------------- Check Put-Call Parity -------------" << std :: endl ;
    auto parity_lhs = vanillaCallPricer.value()-vanillaPutPricer.value();
    auto parity_rhs = S0 - (domDisc/forDisc)*K;
    auto parity_spread = std::abs(parity_lhs - parity_rhs);
    std :: cout << "Call - Put " << parity_lhs << std::endl;
    std :: cout << "Spot - Discounted strike " << parity_rhs << std::endl;
    std :: cout << "Put-Call Parity Spread " << parity_spread << std::endl;
    std :: cout << " -------------- AON Values -------------" << std :: endl ;
    std :: cout << " Value :" << aonCallPricer . value () << std :: endl ;
    std :: cout << " Delta :" << aonCallPricer . delta () << std :: endl ;
    std :: cout << " Gamma :" << aonCallPricer . gamma () << std :: endl ;
    std :: cout << " Vega :" << aonCallPricer . vega ( tau ) << std :: endl ;
    std :: cout << " Theta :" << aonCallPricer . theta ( tau ) << std :: endl ;
    std :: cout << " Delta Fwd :" << aonCallPricer . deltaForward () << std :: endl ;
    std :: cout << " Gamma Fwd :" << aonCallPricer . gammaForward () << std :: endl ;


}

int main(int, char**) {
    std::cout << "Hello, world!\n";
    // Date Stuff
    QuantLib::Calendar calendar = QuantLib::TARGET();
    QuantLib::Date todaysDate(6, QuantLib::December, 2018); // today and the settlement date are the same
    ql::Settings::instance().evaluationDate() = todaysDate; // I think it should be set
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
    QuantLib::Handle<QuantLib::Quote> spotPrice(boost::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(stock))); // current stock price
    QuantLib::Handle<QuantLib::YieldTermStructure> flatDividendTS(boost::shared_ptr<QuantLib::YieldTermStructure>(
            new QuantLib::FlatForward(
                settlementDate,
                dividendYield,
                dayCounter)));
    // flat TS for the risk-free rate
    QuantLib::Handle<QuantLib::YieldTermStructure> flatTS(boost::shared_ptr<QuantLib::YieldTermStructure>(
            new QuantLib::FlatForward(
                settlementDate,
                riskFreeRate,
                dayCounter)));
    // flat TS for volatility
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
    boost::shared_ptr<QuantLib::PricingEngine> Engine(new QuantLib::AnalyticEuropeanEngine(bsmProcess)); // we need to know how to price a derivative
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

    testSolver1();
    
    return 0;
 
    
}
