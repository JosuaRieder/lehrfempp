#include "outside.h"

void outside() {
  std::cout << "Outside: x = " << lf::base::cv::Get<int>("x") << "\n";
  std::cout << "Outside: y = " << lf::base::cv::Get<std::string>("y") << "\n";
}
