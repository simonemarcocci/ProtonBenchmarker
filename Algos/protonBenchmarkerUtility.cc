#include "protonBenchmarkerUtility.h"

namespace pbutil{

  bool protonBenchmarkerUtility::isInTPC(const art::Ptr< simb::MCParticle >& a){
    
    float sx = a->Vx();
    float sy = a->Vy();
    float sz = a->Vz();
    //float ex = a->EndX();
    //float ey = a->EndY();
    //float ez = a->EndZ();

    if ( (sx > 0 && sx < 256 && sy > -111.5 && sy < 111.5 && sz > 0 && sz < 1036) 
        /*|| (ex > 0 && ex < 256 && ey > -111.5 && ey < 111.5 && ez > 0 && ez < 1036)*/)
        return true;
    else return false;

  }
  
  bool protonBenchmarkerUtility::isInTPC(const simb::MCParticle & a){
    
    float sx = a.Vx();
    float sy = a.Vy();
    float sz = a.Vz();
    //float ex = a->EndX();
    //float ey = a->EndY();
    //float ez = a->EndZ();

    if ( (sx > 0 && sx < 256 && sy > -111.5 && sy < 111.5 && sz > 0 && sz < 1036) 
        /*|| (ex > 0 && ex < 256 && ey > -111.5 && ey < 111.5 && ez > 0 && ez < 1036)*/)
        return true;
    else return false;

  }
  

  std::vector<double> protonBenchmarkerUtility::getMomentumVector(const art::Ptr< simb::MCParticle >& a){

    std::vector<double> aMom;
    aMom = {a->Px(), a->Py(), a->Pz()};

    return aMom;

  }

  std::vector<double> protonBenchmarkerUtility::getMomentumVector(const recob::Track& a){

    std::vector<double> aMom;
    auto const& smv = a.StartMomentumVector();
    aMom = {smv.X(), smv.Y(), smv.Z()};

    return aMom;

  }

  float protonBenchmarkerUtility::getAngle(const std::vector<double> a, const std::vector<double> b, protonBenchmarkerUtility rbutil, std::string proj){

    std::vector<double> aMom;
    std::vector<double> bMom;

    if (proj == "no"){

      aMom = {a.at(0), a.at(1), a.at(2)};
      bMom = {b.at(0), b.at(1), b.at(2)};

    }

    else if ((proj == "xz") | (proj == "zx")){

      aMom = {a.at(0), a.at(2)};
      bMom = {b.at(0), b.at(2)};

    }

    else if ((proj == "xy") | (proj == "yx")){

      aMom = {a.at(0), a.at(1)};
      bMom = {b.at(0), b.at(1)};

    }

    else if ((proj == "yz") | (proj == "zy")){

      aMom = {a.at(1), a.at(2)};
      bMom = {b.at(1), b.at(2)};

    }
    else
      throw std::invalid_argument("Valid arguments are 'no', 'xy', 'xz' and 'yz'");


    std::vector<float> aMomUnit = rbutil.getUnitVector(aMom);
    std::vector<float> bMomUnit = rbutil.getUnitVector(bMom);

    float angle = std::acos(rbutil.getDotProduct(aMomUnit, bMomUnit)) * 180 / 3.14159;

    return angle;
  }

  std::vector<float> protonBenchmarkerUtility::getUnitVector(std::vector<double> a){
/*
    float ax = a.at(0);
    float ay = a.at(1);
    float az = a.at(2);

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

  float protonBenchmarkerUtility::getDotProduct(std::vector<float> a, std::vector<float> b){

    if (a.size() != b.size())
      throw std::invalid_argument("Cannot dot product vectors of different sizes");
  
    float dotProduct = 0;
    for (size_t i = 0; i < a.size(); i++){

      dotProduct += a.at(i)*b.at(i);

    }

    return dotProduct;

  }
 
  std::vector<float> protonBenchmarkerUtility::getHitXZPosition(const recob::Hit& thisHit, protonBenchmarkerUtility rbutil){

    float hitChannel = (float)thisHit.Channel();
    float hitX = rbutil.convertTicksToX(thisHit);
    float hitZ = (hitChannel-4800)*0.3;
      
    std::vector<float> hitPosition = {hitX, hitZ};

    return hitPosition;

  }

  float protonBenchmarkerUtility::convertTicksToX(const recob::Hit& thisHit){

    double tick = thisHit.PeakTime();

    const detinfo::DetectorProperties* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    //geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
   

    geo::PlaneID planeID;

    if (thisHit.Channel() < 2400)
      planeID = geo::PlaneID(0,0,0);
    else if (thisHit.Channel() >= 4800)
      planeID = geo::PlaneID(0,0,2);
    else
      planeID = geo::PlaneID(0,0,1);

    float xPos = detProp->ConvertTicksToX(tick, planeID);


    return xPos;

  }

  bool protonBenchmarkerUtility::isHitNearVertex(std::vector<float> v, std::vector<float> h){

    float maxDistanceFromVertex = 5.0;

    float vx = v.at(0);
    float vz = v.at(1);
    float hx = h.at(0);
    float hz = h.at(1);

    //std::cout << "v: " << vx << "," << vz << "    h: " << hx << "," << hz << std::endl;

    float hvLength = sqrt(std::pow(vx-hx,2) + std::pow(vz-hz,2));

    if (std::abs(hvLength) <  maxDistanceFromVertex)
      return true;
    else 
      return false;
      
  }

}
