#include "recoBenchmarkerUtility.h"

namespace rbutil{

  float recoBenchmarkerUtility::getAngle(const art::Ptr<simb::MCParticle>& a, const art::Ptr<simb::MCParticle>& b, recoBenchmarkerUtility rbutil, std::string proj){

    std::vector<double> aMom;
    std::vector<double> bMom;

    if (proj == "no"){

      aMom = {a->Px(), a->Py(), a->Pz()};
      bMom = {b->Px(), b->Py(), b->Pz()};

    }

    else if ((proj == "xz") | (proj == "zx")){

      aMom = {a->Px(), a->Pz()};
      bMom = {b->Px(), b->Pz()};

    }

    else if ((proj == "xy") | (proj == "yx")){

      aMom = {a->Px(), a->Py()};
      bMom = {b->Px(), b->Py()};

    }

    else if ((proj == "yz") | (proj == "zy")){

      aMom = {a->Py(), a->Pz()};
      bMom = {b->Py(), b->Pz()};

    }
    else
      throw std::invalid_argument("Valid arguments are 'no', 'xy', 'xz' and 'yz'");


    std::vector<float> aMomUnit = rbutil.getUnitVector(aMom);
    std::vector<float> bMomUnit = rbutil.getUnitVector(bMom);

    float angle = std::acos(rbutil.getDotProduct(aMomUnit, bMomUnit)) * 180 / 3.14159;

    return angle;
  }

  std::vector<float> recoBenchmarkerUtility::getUnitVector(std::vector<double> a){
/*
    float ax = a->Px();
    float ay = a->Py();
    float az = a->Pz();

    float aMag = std::sqrt(ax*ax + ay*ay + az*az);

    float axNorm = ax / aMag;
    float ayNorm = ay / aMag;
    float azNorm = az / aMag;
*/

    std::vector<float> aNorm;

    double magnitude = 0;
    for (size_t i = 0; i < a.size(); i++){

      magnitude += (a.at(i) * a.at(i));

    }

    magnitude = std::sqrt(magnitude);

    for (size_t i = 0; i < a.size(); i++){

      aNorm.push_back(a.at(i)/magnitude);

    }

    return aNorm;

  }

  float recoBenchmarkerUtility::getDotProduct(std::vector<float> a, std::vector<float> b){

    if (a.size() != b.size())
      throw std::invalid_argument("Cannot dot product vectors of different sizes");
  
    float dotProduct = 0;
    for (size_t i = 0; i < a.size(); i++){

      dotProduct += a.at(i)*b.at(i);

    }

    return dotProduct;

  }

}
