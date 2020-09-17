#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <vector>
#include <array>

double pi = 3.14159265358979323846;

class Complex{
public:
  double Re, Im;

  Complex(){
    Re = 0; Im = 0;
  };
  Complex(double a, double b){
    Re = a; Im = b;
  };
  double Get_Re(){
    return Re;
  };
  double Get_Im(){
    return Im;
  };
  double Abs(){
    return std::sqrt(Re*Re+Im*Im);
  };
  void Conjugate(){
    Im = -1 * Im;
  };

};

Complex Add(Complex a, Complex b){
  Complex c;
  c.Re = a.Re + b.Re;
  c.Im = a.Im + b.Im;
  return c;  
}

Complex Subtract(Complex a, Complex b){
  Complex c;
  c.Re = a.Re - b.Re;
  c.Im = a.Im - b.Im;
  return c;  
}

Complex Multiply(Complex a, Complex b){
  Complex c;
  c.Re = a.Re * b.Re - a.Im * b.Im;
  c.Im = a.Re * b.Im + a.Im * b.Re;
  return c;
}

Complex Multiply(double a, Complex b){
  Complex c;
  c.Re = b.Re * a;
  c.Im = b.Im * a;
  return c;
}

Complex Conjugate(Complex a){
  Complex c;
  c.Im = a.Im * (-1);
  return c;
}

//////////////////////////////
std::string decimal_to_binary(int N){
  if(N == 0) return "0";
  int bin_digit = std::floor(std::log2(N)) + 1;
  std::string str = "";
  int q;
  for(int i = 0 ; i < bin_digit ; i++){
    int p = std::pow(2, bin_digit - 1 - i);
    q = N / p;
    if(q == 1){
      str += '1';
      N -= p;
    }else{
      str += '0';
    }
  }
  return str;
}

std::string decimal_to_binary(int N , int digit){
  if(N == 0) return "0";
  int bin_digit = std::floor(std::log2(N)) + 1;
  if(bin_digit <= digit){
    std::string str = "";
    int q;
    for(int i = 0 ; i < digit ; i++){
      int p = std::pow(2, digit - 1 - i);
      q = N / p;
      if(q == 1){
	str += '1';
	N -= p;
      }else{
	str += '0';
      }
    }
    return str;
  }else{
    return "error";
  }
}

std::string binary_inversion(int N, int digit){
  std::string str = decimal_to_binary(N, digit);
  std::string str_inv = "";
  for(int i = (int)str.size() - 1; i >= 0 ; i--){
    str_inv += str[i];
  }
  return str_inv;
}

int binary_to_decimal(std::string binary_string){
  int a = 0;
  for(int i = 0 ; i < (int)binary_string.size() ; i++){
    if(binary_string[i] == '1'){
      a += pow(2, (double)binary_string.size() - i - 1);
    }else{
      continue;
    }
  }
  return a;
}


void FFT(std::vector<double> data, std::vector<double> &transform_Re, std::vector<double> &transform_Im){
  if(data.size() != transform_Re.size() || data.size() != transform_Im.size()){
    return;
  }

  int N = data.size(); // 2^n
  int n = std::log2(N);
  std::vector<int> seq(N);
  std::iota(seq.begin(), seq.end(), 0.0);
  std::vector<int> bit_inv(N);
  for(int i = 0 ; i < N ; i++){
    bit_inv[i] = binary_to_decimal(binary_inversion(seq[i] , n));
  }

  std::vector<Complex> Cdata(N); // Data as a complex number
  for(int i = 0 ; i < N ; i++){  // Initialization;
    Cdata[i].Re = data[bit_inv[i]];
    Cdata[i].Im = 0;
  }

  
  std::vector<Complex> aft_btfly(N); // The result of butterfly operation
  for(int i = 1 ; i <= n ; i++){
    //std::cout<<"i= "<<i<<std::endl;
    int size = std::pow(2,i);
    std::vector<Complex> W(size / 2);
    //std::cout << "size = " << W.size() << std::endl;
    for(int l = 0 ; l < (int)W.size() ; l++){
      W[l].Re = std::cos(l * 2 * pi / size);
      W[l].Im = std::sin(l * 2 * pi / size) * (-1);
      //std::cout << W[l].Re << " " << W[l].Im << std::endl;;
    }
    
    for(int k = 0 ; k < (int)N / std::pow(2, i) ; k++){
      for(int j = 0 ; j < (int)std::pow(2, i - 1) ; j++){

	//aft_btfly[size * k + j]          = data[ size * k + j] ] + W[j] * data[ size * k + size/2 + j]];
	//aft_bufly[size * k + size/2 + j] = data[ size * k + j] ] - W[j] * data[ size * k + size/2 + j]];
	
	aft_btfly[size * k + j]          = Add     (Cdata[size * k + j], Multiply(W[j], Cdata[size * k + size/2 + j]));
	aft_btfly[size * k + size/2 + j] = Subtract(Cdata[size * k + j], Multiply(W[j], Cdata[size * k + size/2 + j]));
	
      }//j
    }//k
    for(int ii = 0 ; ii < N ; ii++){
      Cdata[ii] = aft_btfly[ii];
    }
  }//i
  for(int i = 0 ; i < N ; i++){
    transform_Re[i] = Cdata[i].Re;
    transform_Im[i] = Cdata[i].Im;
  }
  return;
}


void InverseFFT(std::vector<double> Re, std::vector<double> Im, std::vector<double> &transform_Re, std::vector<double> &transform_Im){
  int N = Re.size(); // 2^n
  int n = std::log2(N);
  std::vector<int> seq(N);
  std::iota(seq.begin(), seq.end(), 0.0);
  std::vector<int> bit_inv(N);
  for(int i = 0 ; i < N ; i++){
    bit_inv[i] = binary_to_decimal(binary_inversion(seq[i] , n));
  }

  std::vector<Complex> Cdata(N); // Data as a complex number
  for(int i = 0 ; i < N ; i++){  // Initialization;
    Cdata[i].Re = Re[bit_inv[i]];
    Cdata[i].Im = Im[bit_inv[i]] * (-1);
  }

  
  std::vector<Complex> aft_btfly(N); // The result of butterfly operation
  for(int i = 1 ; i <= n ; i++){
    //std::cout<<"i= "<<i<<std::endl;
    int size = std::pow(2,i);
    std::vector<Complex> W(size / 2);
    //std::cout << "size = " << W.size() << std::endl;
    for(int l = 0 ; l < (int)W.size() ; l++){
      W[l].Re = std::cos(l * 2 * pi / size);
      W[l].Im = std::sin(l * 2 * pi / size) * (-1);
      //std::cout << W[l].Re << " " << W[l].Im << std::endl;;
    }
    
    for(int k = 0 ; k < (int)N / std::pow(2, i) ; k++){
      for(int j = 0 ; j < (int)std::pow(2, i - 1) ; j++){

	//aft_btfly[size * k + j]          = data[ size * k + j] ] + W[j] * data[ size * k + size/2 + j]];
	//aft_bufly[size * k + size/2 + j] = data[ size * k + j] ] - W[j] * data[ size * k + size/2 + j]];
	
	aft_btfly[size * k + j]          = Add     (Cdata[size * k + j], Multiply(W[j], Cdata[size * k + size/2 + j]));
	aft_btfly[size * k + size/2 + j] = Subtract(Cdata[size * k + j], Multiply(W[j], Cdata[size * k + size/2 + j]));
	
      }//j
    }//k
    for(int ii = 0 ; ii < N ; ii++){
      Cdata[ii] = aft_btfly[ii];
    }
  }//i
  for(int i = 0 ; i < N ; i++){
    transform_Re[i] = Cdata[i].Re / N;
    transform_Im[i] = Cdata[i].Im * (-1) / N;
  }
  return;  
}
