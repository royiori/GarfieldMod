#pragma once

#include "ComponentFieldMap.hh"

namespace Garfield {

/// Component for importing and interpolating Comsol field maps.

class ComponentComsol : public ComponentFieldMap {
 public:
  /// Default constructor
  ComponentComsol();
  ComponentComsol(std::string mesh, std::string mplist, std::string field);
  /// Destructor
  ~ComponentComsol() {}

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  void ElectricField2(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status);

  void SetGridXYZ(double _nx, double _ny, double _nz) {nx=_nx; ny=_ny; nz=_nz;}


  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;

  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  Medium* GetMedium(const double x, const double y, const double z) override;

  bool Initialise(std::string header = "mesh.mphtxt",
                  std::string mplist = "dielectrics.dat",
                  std::string field = "field.txt");

  // Two methods to add the weighting field --LIUQ
  //   1. generate a seperate weighting field file in COMSOL, and read by SetWeightingField( filename, label ); 
  //   2. the weighting fields & orginal field are all stored in the "field.txt", and will be read in automatically by Initialise without label.
  //      so call SetWeightingFieldLabel( id, label ) to add a label attached to each of them.
  //      It's user's responsibility to make sure the id & label matched with each other.
  //
  bool SetWeightingField(std::string file, std::string label);
  bool SetWeightingFieldLabel(int id, std::string label); 

 protected:
  void UpdatePeriodicity() override { UpdatePeriodicityCommon(); }

  double GetElementVolume(const unsigned int i) override;
  void GetAspectRatio(const unsigned int i, double& dmin,
                      double& dmax) override;

 
  bool mapField; //LIUQ
  double xmin, xmax, xstep;
  double ymin, ymax, ystep;
  double zmin, zmax, zstep;
  double nx, ny, nz;

  std::vector<double>xlist;
  std::vector<double>ylist;
  std::vector<double>zlist;
  std::vector<double>Exlist;
  std::vector<double>Eylist;
  std::vector<double>Ezlist;
  std::vector<double>vlist;
  std::vector<Medium*>medlist;
  std::vector<int>statlist;

  struct nodeCmp {
    bool operator()(const ComponentFieldMap::Node& lhs,
                    const ComponentFieldMap::Node& rhs) const {
      double dx = round(lhs.x * 1e6) - round(rhs.x * 1e6);
      double dy = round(lhs.y * 1e6) - round(rhs.y * 1e6);
      double dz = round(lhs.z * 1e6) - round(rhs.z * 1e6);
      return dx < 0 || (dx == 0 && (dy < 0 || (dy == 0 && dz < 0)));
    }
  };
};
}
