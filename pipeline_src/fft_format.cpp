#include "fft_format.h"

std::string formatDouble7(double value) {
  std::ostringstream oss;

  if (value >= 100000 || value <= -10000) {
    return " ERROR ";
  }

  if (std::abs(value) >= 1000) {
    oss << std::fixed << std::setprecision(2);
  } else if (std::abs(value) >= 100) {
    oss << std::fixed << std::setprecision(3);
  } else if (std::abs(value) >= 10) {
    oss << std::fixed << std::setprecision(4);
  } else {
    oss << std::fixed << std::setprecision(5);
  }

  oss << value;

  std::string result = oss.str();

  if (result.length() > 7) {
    result = result.substr(0, 7);
  } else if (result.length() < 7) {
    result.insert(0, 7 - result.length(), ' ');
  }

  return result;
}