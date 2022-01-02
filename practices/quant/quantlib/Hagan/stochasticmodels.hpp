/*
 * stochasticmodels.hpp
 *
 *  Created on: Jan 2, 2022
 *      Author: jaeyong
 */

#ifndef STOCHASTICMODELS_HPP_
#define STOCHASTICMODELS_HPP_

#include <iostream>
 
using namespace std;

// Base class
class Shape {
   public:
      void setWidth(int w) {
         width = w;
      }
      void setHeight(int h) {
         height = h;
      }
      template <typename F> void test(F arg0){
        std::cout << "test function in the Shape class must be implmented in children classes!" << std::endl;
        exit(-1);
      }
      
   protected:
      int width;
      int height;
};

    //   struct drift{
    //     size_t lhs_sv;
    //     double coeff;
    //     size_t rhs_sv;
    //     double order;
    //   };    

// Derived class
class Rectangle: public Shape {
   public:
      int getArea() { 
         return (width * height); 
      };
      struct drift;
      void test(drift arg);
};

struct Rectangle::drift{
    size_t lhs_sv;
    double coeff;
    size_t rhs_sv;
    double order;
    };    

#endif /* STOCHASTICMODELS_HPP_ */
