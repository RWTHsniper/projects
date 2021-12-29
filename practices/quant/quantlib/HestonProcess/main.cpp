#include <iostream>
#include <ql/quantlib.hpp>


namespace ql = QuantLib;

void testingStochasticProcesses3 (){
    ql::Date refDate = ql::Date (27 , ql::Sep ,2009);
    ql::Rate riskFreeRate =0.0321;
    ql::Rate dividendRate =0.0128;
    ql::Real spot =52.0;
    ql::Rate vol =0.2144;
    ql::Calendar cal = ql::TARGET ();
    ql::DayCounter dc= ql::Actual365Fixed ();
    boost :: shared_ptr < ql::YieldTermStructure > rdStruct ( new ql::FlatForward ( refDate , riskFreeRate ,dc ));
    boost :: shared_ptr < ql::YieldTermStructure > rqStruct ( new ql::FlatForward ( refDate , dividendRate ,dc ));
    ql::Handle < ql::YieldTermStructure > rdHandle ( rdStruct );
    ql::Handle < ql::YieldTermStructure > rqHandle ( rqStruct );
    boost :: shared_ptr < ql::SimpleQuote > spotQuote (new ql::SimpleQuote ( spot ));
    ql::Handle <ql::Quote > spotHandle ( spotQuote );
    boost :: shared_ptr < ql::BlackVolTermStructure > volQuote (new ql::BlackConstantVol ( refDate , cal , vol , dc ));
    ql::Handle < ql::BlackVolTermStructure > volHandle ( volQuote );
    ql::Real v0 =0.12 , kappa =1.2 , theta =0.08 , sigma =0.05 , rho = -0.6;
    boost :: shared_ptr < ql::HestonProcess > hestonProcess ( new ql::HestonProcess ( rdHandle , rqHandle , spotHandle ,v0 ,
    kappa ,theta ,sigma ,rho , ql::HestonProcess :: PartialTruncation ));
    ql::BigInteger seed =12324;
    ql::MersenneTwisterUniformRng unifMt ( seed );
    ql::BoxMullerGaussianRng < ql::MersenneTwisterUniformRng > bmGauss ( unifMt );
    ql::Time dt =0.10 ,t =0.0;
    ql::Array dw (2) ,x (2);
    // x is the 2- dimensional process
    x [0]= spotQuote -> value ();
    x [1]= v0;
    ql::Size numVals =10;
    for ( ql::Size j=1;j <= numVals ;++j){
        dw [0]= bmGauss . next (). value ;
        dw [1]= bmGauss . next (). value ;
        x= hestonProcess -> evolve (t,x,dt ,dw );
        std :: cout << " Time : " << t+dt << ", S_t : " << x[0] << ", V_t : " << x [1] << std :: endl ;
        t+= dt;
    }
}

int main(int, char**) {
    std::cout << "Hello, world!\n";
    testingStochasticProcesses3();
}