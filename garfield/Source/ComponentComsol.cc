// Copied and modified ComponentAnsys123.cc
// >>Modified by Qian LIU << 2019.4.15
// >>Add the following features:
// >>
// >>    1. add weighting filed methods
// >>       Two methods to add the weighting field
// >>          1. generate a seperate weighting field file in COMSOL, and read by SetWeightingField( filename, label ); 
// >>          2. the weighting fields & orginal field are all stored in the "field.txt", and will be read in automatically by Initialise without label.
// >>             so call SetWeightingFieldLabel( id, label ) to add a label attached to each of them.
// >>             It's user's responsibility to make sure the id & label matched with each other.
// >>
// >>    2. the comsol file formats are updated to COMSOL 5.3a
// >>
// >>    3. protection for Win/Mac comsol generated files
// >>
// >>    4. protection for the end of file reached
// >>
// >>
// >>
// >>

#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "ComponentComsol.hh"

namespace Garfield {

ComponentComsol::ComponentComsol() : ComponentFieldMap() {
  mapField = false;
  xmin = 0, ymin = 0, zmin = 0;
  xmax = 0, ymax = 0, zmax = 0;
  nx =  40, ny =  40, nz =  40;
  m_className = "ComponentComsol";
}

ComponentComsol::ComponentComsol(std::string mesh, std::string mplist,
                                 std::string field)
    : ComponentFieldMap() {
  mapField = false;
  xmin = 0, ymin = 0, zmin = 0; 
  xmax = 0, ymax = 0, zmax = 0;
  nx =  40, ny =  40, nz =  40;
  m_className = "ComponentComsol";
  Initialise(mesh, mplist, field);
}

bool ends_with(std::string ss, std::string tt) {
  if(ss.length()==0) return 0;
  if(ss.at(ss.length()-1)=='\r') ss.resize(ss.size()-1); //remove the '\r'
  return ss.size() >= tt.size() && ss.substr(ss.size() - tt.size(), tt.size()) == tt;
}

int readInt(std::string s) {
  std::istringstream iss(s);
  int ret;
  iss >> ret;
  return ret;
}

bool ComponentComsol::Initialise(std::string mesh, std::string mplist,
                                 std::string field) {
  m_ready = false;
  m_warning = false;
  m_nWarnings = 0;
  double unit = 100.0;  // m
  std::string line;

  int LIUQ = 1;
  //-------------------------
  //1. Open the materials file.
  materials.clear();
  std::ifstream fmplist;
  fmplist.open(mplist.c_str(), std::ios::in);
  if (fmplist.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open result file " << mplist
              << " for reading.\n";
    return false;
  }
  if(LIUQ) std::cout << "--> Reading "<<mplist<<"."<<std::endl;

  //1.1 read number of materials
  fmplist >> m_nMaterials;

  for (unsigned int i = 0; i < m_nMaterials; ++i) {
    Material newMaterial;
    newMaterial.driftmedium = true;
    newMaterial.medium = nullptr;
    newMaterial.ohm = -1;
    fmplist >> newMaterial.eps;
    if(fabs(newMaterial.eps-1)>1E-3) newMaterial.driftmedium = false;
    materials.push_back(newMaterial);
  }
  if(LIUQ) std::cout << ".\t"<<m_nMaterials<<" materials are read."<<std::endl;
  {
    // add default material
    Material newMaterial;
    newMaterial.driftmedium = false;
    newMaterial.medium = nullptr;
    newMaterial.eps = newMaterial.ohm = -1;
    materials.push_back(newMaterial);
    m_nMaterials++;
  }

  //1.2 read domain & material map
  std::map<int, int> domain2material;
  int d2msize;
  fmplist >> d2msize;
  for (int i = 0; i < d2msize; ++i) {
    int domain;
    fmplist >> domain;
    fmplist >> domain2material[domain];
  }
  if(LIUQ) std::cout << ".\t"<<d2msize<<" domains are read."<<std::endl;
  fmplist.close();


  //-------------------------
  //2. Open the node file
  nodes.clear();
  std::ifstream fmesh;
  fmesh.open(mesh.c_str(), std::ios::in);
  if (fmesh.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open nodes file " << mesh << " for reading.\n";
    return false;
  }
  if(LIUQ) std::cout << "--> Reading "<<mesh<<"."<<std::endl;

  //------
  //2.1 find the filehead of the mesh points
  do {
    std::getline(fmesh, line); 
  } while (!ends_with(line, "# number of mesh points") && fmesh.peek() != EOF);

  if(fmesh.peek()==EOF) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    The " << mesh << " format is not correct, check the COMSOL version.\n";
    return false;
  }
 


  nNodes = readInt(line);
  std::cout << ".\tGoing to read " << nNodes << " nodes from file " << mesh << ".\n";

  //------
  //2.2 read mesh points (X,Y,Z) and store to nodes
  do {
    std::getline(fmesh, line);
  } while (!ends_with(line, "# Mesh point coordinates") && fmesh.peek() != EOF);

  if(fmesh.peek()==EOF) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    The " << mesh << " format is not correct, check the COMSOL version.\n";
    return false;
  }


  xmin = 1e100, ymin = 1e100, zmin = 1e100;
  xmax = -1e100, ymax = -1e100, zmax = -1e100;
  for (int i = 0; i < nNodes; ++i) {
    Node newNode;
    fmesh >> newNode.x >> newNode.y >> newNode.z;
    newNode.x *= unit;
    newNode.y *= unit;
    newNode.z *= unit;
    nodes.push_back(newNode);
    xmin = std::min(xmin, newNode.x);
    xmax = std::max(xmax, newNode.x);
    ymin = std::min(ymin, newNode.y);
    ymax = std::max(ymax, newNode.y);
    zmin = std::min(zmin, newNode.z);
    zmax = std::max(zmax, newNode.z);
  }
  if(LIUQ) std::cout << ".\t" << nodes.size() << " nodes are read."<<std::endl;


  //------
  //2.3 find tet2 meshes
  do {
    std::getline(fmesh, line);
  } while (!ends_with(line, "4 tet2 # type name")  && fmesh.peek() != EOF);

  if(fmesh.peek()==EOF) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    The " << mesh << " format is not correct, check the COMSOL version.\n";
    return false;
  }


  //------
  //2.4 read elements (10 nodes(mesh points) form an element)
  do {
    std::getline(fmesh, line);
  } while (!ends_with(line, "# number of elements") && fmesh.peek() != EOF);

  if(fmesh.peek()==EOF) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    The " << mesh << " format is not correct, check the COMSOL version.\n";
    return false;
  }

  nElements = readInt(line);
  elements.clear();
  std::cout << ".\tGoing to read " << nElements << " elements from file " << mesh
            << ".\n";
  std::getline(fmesh, line);
  // elements 6 & 7 are swapped due to differences in COMSOL and ANSYS
  // representation
  int perm[10] = {0, 1, 2, 3, 4, 5, 7, 6, 8, 9};
  for (int i = 0; i < nElements; ++i) {
    Element newElement;
    newElement.degenerate = false;
    for (int j = 0; j < 10; ++j) {
      fmesh >> newElement.emap[perm[j]];
    }
    elements.push_back(newElement);
  }
  if(LIUQ) std::cout << ".\t" << elements.size() << " elements are read."<<std::endl;


  //------
  // 2.5 read the domain number for each element (domain number is the same as in COMSOL material table)
  //     then set the material to each element according to the domain-material table
  do {
    std::getline(fmesh, line);
  } while (!ends_with(line,"# Geometric entity indices") && fmesh.peek() != EOF);

  if(fmesh.peek()==EOF) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    The " << mesh << " format is not correct, check the COMSOL version.\n";
    return false;
  }

  for (int i = 0; i < nElements; ++i) {
    int domain;
    fmesh >> domain;
    elements[i].matmap = domain2material.count(domain) ? domain2material[domain]
                                                       : m_nMaterials - 1;
  }
  fmesh.close();

  std::map<Node, std::vector<int>, nodeCmp> nodeIdx;
  for (int i = 0; i < nNodes; ++i) {
    nodeIdx[nodes[i]].push_back(i);
  }
  std::cout << ".\tMapping elements to " << nodeIdx.size() << " domain nodes."<<std::endl;





  //-------------------------
  //3. Open the field file
  std::ifstream ffield;
  ffield.open(field.c_str(), std::ios::in);
  if (ffield.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open field potentials file " << field
              << " for reading.\n";
    return false;
  }
  if(LIUQ) std::cout << "--> Reading "<<field<<"."<<std::endl;

  //------
  //3.1 search to the file head
  do {
    std::getline(ffield, line);
  } while (line.substr(0, 81) !=
           "% X                       Y                        Z               "
           "         V (V)");

  //3.2 for field: X, Y, Z, V
  //    for weighting field: X, Y, Z, V0, V1, V2, V3.... nWeightingFields is the total number of voltages
  //{
  //  std::istringstream sline(line);
  //  std::string token;
  //  sline >> token;  // %
  //  sline >> token;  // x
  //  sline >> token;  // y
  //  sline >> token;  // z
  //  sline >> token;  // V
  //  sline >> token;  // (V)
  //  while (sline >> token) {
  //    std::cout << m_className << "::Initialise:\n";
  //    std::cout << "    Reading data for weighting field " << token << ".\n";
  //    nWeightingFields++;
  //    wfields.push_back(token);
  //    wfieldsOk.push_back(true);
  //    sline >> token;  // (V)
  //  }
  //}
  nWeightingFields = -1;
  unsigned long loc=0;
  while((loc=line.find( "V (V) @ ", loc ))!=std::string::npos) 
  {
     nWeightingFields++;
     loc++;
  }
  if(nWeightingFields==-1) nWeightingFields = 0;
  wfields.resize(nWeightingFields);
  wfieldsOk.resize(nWeightingFields);
  if(LIUQ) std::cout<<".\tn="<<nWeightingFields<<" weighting field(s) will be read."<<std::endl;

  //3.3 read the field, and store to the nodes
  for (int i = 0; i < nNodes; ++i) {
    Node tmp;
    ffield >> tmp.x >> tmp.y >> tmp.z >> tmp.v;
    tmp.x *= unit;
    tmp.y *= unit;
    tmp.z *= unit;
    for (int j = 0; j < nWeightingFields; ++j) {
      double w;
      ffield >> w;
      tmp.w.push_back(w);
    }

    // compare the field with the nodes
    int closest = -1;
    double closestDist = 1;
    const unsigned int nIdx = nodeIdx[tmp].size();
    // for (int j : nodeIdx[tmp]) {
    for (unsigned int k = 0; k < nIdx; ++k) {
      int j = nodeIdx[tmp][k];
      double dist = (tmp.x - nodes[j].x) * (tmp.x - nodes[j].x) +
                    (tmp.y - nodes[j].y) * (tmp.y - nodes[j].y) +
                    (tmp.z - nodes[j].z) * (tmp.z - nodes[j].z);
      if (dist < closestDist) {
        closestDist = dist;
        closest = j;
      }
      
    }
    if (closest == -1) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Could not match the node from field potentials file: "
                << tmp.x << " " << tmp.y << " " << tmp.z << "\n.";
      return false;
    }

    nodes[closest].v = tmp.v;
    nodes[closest].w = tmp.w;
  }

  m_ready = true;

  //for (int i = 0; i < nNodes; ++i) {
  //  double ex, ey, ez, v;
  //  Medium* m;
  //  int status;
  //  ElectricField(nodes[i].x, nodes[i].y, nodes[i].z, ex, ey, ez, v, m, status);
  //  std::cout << "Field at " << nodes[i].x << " " << nodes[i].y << " " << nodes[i].z << ": " << ex << " " << ey << " " << ez << " " << v << "\n";
  //}

  // Establish the ranges.
  SetRange();
  UpdatePeriodicity();
  return true;
}

bool ComponentComsol::SetWeightingFieldLabel(int id, std::string label) {
  int iw = nWeightingFields;
  if(m_debug) {
      std::cout << m_className << "::SetWeightingFieldLabel:\n";
      std::cout << "    "<<iw<<" weighting fields are initialized.\n";
  }

  for (int i = nWeightingFields; i--;) {
    if (wfields[i] == label) {
      iw = i;
      break;
    }
  }

  if (id >= nWeightingFields) {
    std::cout << m_className << "::SetWeightingFieldLabel:\n";
    std::cout << "    Only "<<nWeightingFields<<" weighting field(s) stored, can't set id="<<id<<".\n";
    return false;
  }
  if(iw!=nWeightingFields) {
    std::cout << m_className << "::SetWeightingFieldLabel:\n";
    std::cout << "    Existing weighting field " << label << ".\n";
    return false;
  }

  wfields[id] = label;
  wfieldsOk[id] = true;
  return true;
}

bool ComponentComsol::SetWeightingField(std::string field, std::string label) {
  double unit = 100.0;  // m;

  if (!m_ready) {
    std::cerr << m_className << "::SetWeightingField:\n";
    std::cerr << "    No valid field map is present.\n";
    std::cerr << "    Weighting field cannot be added.\n";
    return false;
  }

  // Open the voltage list.
  std::ifstream ffield;
  ffield.open(field.c_str(), std::ios::in);
  if (ffield.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open field potentials file " << field
              << " for reading.\n";
    return false;
  }

  // Check if a weighting field with the same label alm_ready exists.
  int iw = nWeightingFields;
  for (int i = nWeightingFields; i--;) {
    if (wfields[i] == label) {
      iw = i;
      break;
    }
  }
  if (iw == nWeightingFields) {
    ++nWeightingFields;
    wfields.resize(nWeightingFields);
    wfieldsOk.resize(nWeightingFields);
    for (int j = 0; j < nNodes; ++j) {
      nodes[j].w.resize(nWeightingFields);
    }
  } else {
    std::cout << m_className << "::SetWeightingField:\n";
    std::cout << "    Replacing existing weighting field " << label << ".\n";
  }
  wfields[iw] = label;
  wfieldsOk[iw] = false;
  std::map<Node, std::vector<int>, nodeCmp> nodeIdx;
  for (int i = 0; i < nNodes; ++i) {
    nodeIdx[nodes[i]].push_back(i);
  }
  std::cout << "Map size: " << nodeIdx.size() << std::endl;

  std::string line;
  do {
    std::getline(ffield, line);
  } while (line !=
           "% x                       y                        z               "
           "         V (V)");
  for (int i = 0; i < nNodes; ++i) {
    Node tmp;
    ffield >> tmp.x >> tmp.y >> tmp.z >> tmp.v;
    tmp.x *= unit;
    tmp.y *= unit;
    tmp.z *= unit;
    int closest = -1;
    double closestDist = 1;
    const unsigned int nIdx = nodeIdx[tmp].size();
    // for (int j : nodeIdx[tmp]) {
    for (unsigned int k = 0; k < nIdx; ++k) {
      int j = nodeIdx[tmp][k];
      double dist = (tmp.x - nodes[j].x) * (tmp.x - nodes[j].x) +
                    (tmp.y - nodes[j].y) * (tmp.y - nodes[j].y) +
                    (tmp.z - nodes[j].z) * (tmp.z - nodes[j].z);
      if (dist < closestDist) {
        closestDist = dist;
        closest = j;
      }
    }
    if (closest == -1) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Could not match the node from field potentials file: "
                << tmp.x << " " << tmp.y << " " << tmp.z << "\n.";
      return false;
    }
    nodes[closest].w[iw] = tmp.v;
  }

  return true;
}

void ComponentComsol::ElectricField(const double x, const double y,
                                    const double z, double& ex, double& ey,
                                    double& ez, Medium*& m, int& status) {
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, m, status);
}



void ComponentComsol::ElectricField2(const double xin, const double yin,
                                    const double zin, double& ex, double& ey,
                                    double& ez, double& volt, Medium*& m,
                                    int& status) {
   if(!mapField) { // map field according to user defined grids
     if(xmax == 0 && ymax == 0 && zmax == 0 && xmin == 0 && ymin == 0 && zmin == 0) return;

     xstep = (xmax-xmin)/nx;
     ystep = (ymax-ymin)/ny;
     zstep = (zmax-zmin)/nz;

     double _ex, _ey, _ez, _vol;
     Medium* _m;
     int _stat;

     for(int i=0; i<nx; i++) {
        for(int j=0; j<ny; j++) {
           for(int k=0; k<nz; k++) {
                double xt = xmin + i * xstep + xstep/2;
                double yt = ymin + j * ystep + ystep/2;
                double zt = zmin + k * zstep + zstep/2;
                ElectricField2(xt, yt, zt, _ex, _ey, _ez, _vol, _m, _stat);
                xlist.push_back(xt);
                ylist.push_back(yt);
                zlist.push_back(zt);
                Exlist.push_back(_ex);
                Eylist.push_back(_ey);
                Ezlist.push_back(_ez);
                vlist.push_back(_vol);
                medlist.push_back(_m);
                statlist.push_back(_stat);
           }
        }
     }
     std::cout<<" x is set from ["<<xmin<<", "<<xmax<<"] for "<<nx<<" with step "<<xstep<<std::endl;
     std::cout<<" y is set from ["<<ymin<<", "<<ymax<<"] for "<<ny<<" with step "<<ystep<<std::endl;
     std::cout<<" z is set from ["<<zmin<<", "<<zmax<<"] for "<<nz<<" with step "<<zstep<<std::endl;

     mapField=true;
   } 
   else { // read value according to the mapped field
     // Copy the coordinates
     double x = xin, y = yin, z = zin;

     // Map the coordinates onto field map coordinates
     bool xmirr, ymirr, zmirr;
     double rcoordinate, rotation;
     MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

     // Initial values
     ex = ey = ez = volt = 0.;
     status = 0;
     m = NULL;

     int ix = (x - xmin - xstep/2) / xstep;
     int iy = (y - ymin - ystep/2) / ystep;
     int iz = (z - zmin - zstep/2) / zstep;
     int id = ix*ny*nz + iy*nz + iz;

     //check if it's correct or not
     double _x = xlist[id];
     double _y = ylist[id];
     double _z = zlist[id];
     if(fabs(_x-x)>xstep || fabs(_y-y)>ystep || fabs(_z-z)>zstep) {
         std::cout<<"x is out of range for "<<x<<", in list is "<<_x<<" for "<<ix<<" "<<id<<std::endl;
         std::cout<<"y is out of range for "<<y<<", in list is "<<_y<<" for "<<iy<<" "<<id<<std::endl;
         std::cout<<"z is out of range for "<<z<<", in list is "<<_z<<" for "<<iz<<" "<<id<<std::endl;
     }

     ex = Exlist[id];
     ey = Eylist[id];
     ez = Ezlist[id];
     volt = vlist[id];
     m = medlist[id];
     status = statlist[id];

     if(status<0) return;
     UnmapFields(ex, ey, ez, x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);
     
   }
}





void ComponentComsol::ElectricField(const double xin, const double yin,
                                    const double zin, double& ex, double& ey,
                                    double& ez, double& volt, Medium*& m,
                                    int& status) {
  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool xmirr, ymirr, zmirr;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  // Initial values
  ex = ey = ez = volt = 0.;
  status = 0;
  m = NULL;

  // Do not proceed if not properly initialised.
  if (!m_ready) {
    status = -10;
    PrintNotReady("ElectricField");
    return;
  }

  if (m_warning) PrintWarning("ElectricField");

  // Find the element that contains this point
  double t1, t2, t3, t4, jac[4][4], det;
  const int imap = FindElement13(x, y, z, t1, t2, t3, t4, jac, det);
  if (imap < 0) {
    if (m_debug) {
      std::cout << m_className << "::ElectricField:\n";
      std::cout << "    Point (" << x << ", " << y << ", " << z
                << " not in the mesh.\n";
    }
    status = -6;
    return;
  }

  const Element& element = elements[imap];
  if (m_debug) {
    PrintElement("ElectricField", x, y, z, t1, t2, t3, t4, element, 10);
  }
  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];
  const Node& n4 = nodes[element.emap[4]];
  const Node& n5 = nodes[element.emap[5]];
  const Node& n6 = nodes[element.emap[6]];
  const Node& n7 = nodes[element.emap[7]];
  const Node& n8 = nodes[element.emap[8]];
  const Node& n9 = nodes[element.emap[9]];
  // Tetrahedral field
  volt = n0.v * t1 * (2 * t1 - 1) + n1.v * t2 * (2 * t2 - 1) +
         n2.v * t3 * (2 * t3 - 1) + n3.v * t4 * (2 * t4 - 1) +
         4 * n4.v * t1 * t2 + 4 * n5.v * t1 * t3 + 4 * n6.v * t1 * t4 +
         4 * n7.v * t2 * t3 + 4 * n8.v * t2 * t4 + 4 * n9.v * t3 * t4;
  ex = -(n0.v * (4 * t1 - 1) * jac[0][1] + n1.v * (4 * t2 - 1) * jac[1][1] +
         n2.v * (4 * t3 - 1) * jac[2][1] + n3.v * (4 * t4 - 1) * jac[3][1] +
         n4.v * (4 * t2 * jac[0][1] + 4 * t1 * jac[1][1]) +
         n5.v * (4 * t3 * jac[0][1] + 4 * t1 * jac[2][1]) +
         n6.v * (4 * t4 * jac[0][1] + 4 * t1 * jac[3][1]) +
         n7.v * (4 * t3 * jac[1][1] + 4 * t2 * jac[2][1]) +
         n8.v * (4 * t4 * jac[1][1] + 4 * t2 * jac[3][1]) +
         n9.v * (4 * t4 * jac[2][1] + 4 * t3 * jac[3][1])) /
       det;
  ey = -(n0.v * (4 * t1 - 1) * jac[0][2] + n1.v * (4 * t2 - 1) * jac[1][2] +
         n2.v * (4 * t3 - 1) * jac[2][2] + n3.v * (4 * t4 - 1) * jac[3][2] +
         n4.v * (4 * t2 * jac[0][2] + 4 * t1 * jac[1][2]) +
         n5.v * (4 * t3 * jac[0][2] + 4 * t1 * jac[2][2]) +
         n6.v * (4 * t4 * jac[0][2] + 4 * t1 * jac[3][2]) +
         n7.v * (4 * t3 * jac[1][2] + 4 * t2 * jac[2][2]) +
         n8.v * (4 * t4 * jac[1][2] + 4 * t2 * jac[3][2]) +
         n9.v * (4 * t4 * jac[2][2] + 4 * t3 * jac[3][2])) /
       det;
  ez = -(n0.v * (4 * t1 - 1) * jac[0][3] + n1.v * (4 * t2 - 1) * jac[1][3] +
         n2.v * (4 * t3 - 1) * jac[2][3] + n3.v * (4 * t4 - 1) * jac[3][3] +
         n4.v * (4 * t2 * jac[0][3] + 4 * t1 * jac[1][3]) +
         n5.v * (4 * t3 * jac[0][3] + 4 * t1 * jac[2][3]) +
         n6.v * (4 * t4 * jac[0][3] + 4 * t1 * jac[3][3]) +
         n7.v * (4 * t3 * jac[1][3] + 4 * t2 * jac[2][3]) +
         n8.v * (4 * t4 * jac[1][3] + 4 * t2 * jac[3][3]) +
         n9.v * (4 * t4 * jac[2][3] + 4 * t3 * jac[3][3])) /
       det;

  // Transform field to global coordinates
  UnmapFields(ex, ey, ez, x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);
  //  std::cout << "ef @(" << xin << ", " << yin << ", " << zin << ") = " <<
  // volt << "\n";

  // Drift medium?
  if (m_debug) {
    std::cout << m_className << "::ElectricField:\n";
    std::cout << "    Material " << element.matmap << ", drift flag "
              << materials[element.matmap].driftmedium << "\n";
  }
  m = materials[element.matmap].medium;
  status = -5;
  if (materials[element.matmap].driftmedium) {
    if (m && m->IsDriftable()) status = 0;
  }
}

void ComponentComsol::WeightingField(const double xin, const double yin,
                                     const double zin, double& wx, double& wy,
                                     double& wz, const std::string& label) {
  // Initial values
  wx = wy = wz = 0;

  // Do not proceed if not properly initialised.
  if (!m_ready) return;

  // Look for the label.
  int iw = 0;
  bool found = false;
  for (int i = nWeightingFields; i--;) {
    if (wfields[i] == label) {
      iw = i;
      found = true;
      break;
    }
  }

  // Do not proceed if the requested weighting field does not exist.
  if (!found) return;
  // Check if the weighting field is properly initialised.
  if (!wfieldsOk[iw]) return;

  // Copy the coordinates.
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool xmirr, ymirr, zmirr;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  if (m_warning) PrintWarning("WeightingField");

  // Find the element that contains this point.
  double t1, t2, t3, t4, jac[4][4], det;
  const int imap = FindElement13(x, y, z, t1, t2, t3, t4, jac, det);
  // Check if the point is in the mesh.
  if (imap < 0) return;

  const Element& element = elements[imap];
  if (m_debug) {
    PrintElement("WeightingField", x, y, z, t1, t2, t3, t4, element, 10, iw);
  }
  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];
  const Node& n4 = nodes[element.emap[4]];
  const Node& n5 = nodes[element.emap[5]];
  const Node& n6 = nodes[element.emap[6]];
  const Node& n7 = nodes[element.emap[7]];
  const Node& n8 = nodes[element.emap[8]];
  const Node& n9 = nodes[element.emap[9]];
  // Tetrahedral field
  wx = -(n0.w[iw] * (4 * t1 - 1) * jac[0][1] +
         n1.w[iw] * (4 * t2 - 1) * jac[1][1] +
         n2.w[iw] * (4 * t3 - 1) * jac[2][1] +
         n3.w[iw] * (4 * t4 - 1) * jac[3][1] +
         n4.w[iw] * (4 * t2 * jac[0][1] + 4 * t1 * jac[1][1]) +
         n5.w[iw] * (4 * t3 * jac[0][1] + 4 * t1 * jac[2][1]) +
         n6.w[iw] * (4 * t4 * jac[0][1] + 4 * t1 * jac[3][1]) +
         n7.w[iw] * (4 * t3 * jac[1][1] + 4 * t2 * jac[2][1]) +
         n8.w[iw] * (4 * t4 * jac[1][1] + 4 * t2 * jac[3][1]) +
         n9.w[iw] * (4 * t4 * jac[2][1] + 4 * t3 * jac[3][1])) /
       det;

  wy = -(n0.w[iw] * (4 * t1 - 1) * jac[0][2] +
         n1.w[iw] * (4 * t2 - 1) * jac[1][2] +
         n2.w[iw] * (4 * t3 - 1) * jac[2][2] +
         n3.w[iw] * (4 * t4 - 1) * jac[3][2] +
         n4.w[iw] * (4 * t2 * jac[0][2] + 4 * t1 * jac[1][2]) +
         n5.w[iw] * (4 * t3 * jac[0][2] + 4 * t1 * jac[2][2]) +
         n6.w[iw] * (4 * t4 * jac[0][2] + 4 * t1 * jac[3][2]) +
         n7.w[iw] * (4 * t3 * jac[1][2] + 4 * t2 * jac[2][2]) +
         n8.w[iw] * (4 * t4 * jac[1][2] + 4 * t2 * jac[3][2]) +
         n9.w[iw] * (4 * t4 * jac[2][2] + 4 * t3 * jac[3][2])) /
       det;

  wz = -(n0.w[iw] * (4 * t1 - 1) * jac[0][3] +
         n1.w[iw] * (4 * t2 - 1) * jac[1][3] +
         n2.w[iw] * (4 * t3 - 1) * jac[2][3] +
         n3.w[iw] * (4 * t4 - 1) * jac[3][3] +
         n4.w[iw] * (4 * t2 * jac[0][3] + 4 * t1 * jac[1][3]) +
         n5.w[iw] * (4 * t3 * jac[0][3] + 4 * t1 * jac[2][3]) +
         n6.w[iw] * (4 * t4 * jac[0][3] + 4 * t1 * jac[3][3]) +
         n7.w[iw] * (4 * t3 * jac[1][3] + 4 * t2 * jac[2][3]) +
         n8.w[iw] * (4 * t4 * jac[1][3] + 4 * t2 * jac[3][3]) +
         n9.w[iw] * (4 * t4 * jac[2][3] + 4 * t3 * jac[3][3])) /
       det;

  // Transform field to global coordinates
  UnmapFields(wx, wy, wz, x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);
}

double ComponentComsol::WeightingPotential(const double xin, const double yin,
                                           const double zin,
                                           const std::string& label) {
  // Do not proceed if not properly initialised.
  if (!m_ready) return 0.;

  // Look for the label.
  int iw = 0;
  bool found = false;
  for (int i = nWeightingFields; i--;) {
    if (wfields[i] == label) {
      iw = i;
      found = true;
      break;
    }
  }

  // Do not proceed if the requested weighting field does not exist.
  if (!found) return 0.;
  // Check if the weighting field is properly initialised.
  if (!wfieldsOk[iw]) return 0.;

  // Copy the coordinates.
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates.
  bool xmirr, ymirr, zmirr;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  if (m_warning) PrintWarning("WeightingPotential");

  // Find the element that contains this point.
  double t1, t2, t3, t4, jac[4][4], det;
  const int imap = FindElement13(x, y, z, t1, t2, t3, t4, jac, det);
  if (imap < 0) return 0.;

  const Element& element = elements[imap];
  if (m_debug) {
    PrintElement("WeightingPotential", x, y, z, t1, t2, t3, t4, element, 10,
                 iw);
  }
  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];
  const Node& n4 = nodes[element.emap[4]];
  const Node& n5 = nodes[element.emap[5]];
  const Node& n6 = nodes[element.emap[6]];
  const Node& n7 = nodes[element.emap[7]];
  const Node& n8 = nodes[element.emap[8]];
  const Node& n9 = nodes[element.emap[9]];
  // Tetrahedral field
  return n0.w[iw] * t1 * (2 * t1 - 1) + n1.w[iw] * t2 * (2 * t2 - 1) +
         n2.w[iw] * t3 * (2 * t3 - 1) + n3.w[iw] * t4 * (2 * t4 - 1) +
         4 * n4.w[iw] * t1 * t2 + 4 * n5.w[iw] * t1 * t3 +
         4 * n6.w[iw] * t1 * t4 + 4 * n7.w[iw] * t2 * t3 +
         4 * n8.w[iw] * t2 * t4 + 4 * n9.w[iw] * t3 * t4;
}

Medium* ComponentComsol::GetMedium(const double xin, const double yin,
                                   const double zin) {
  // Copy the coordinates
  double x = xin, y = yin, z = zin;

  // Map the coordinates onto field map coordinates
  bool xmirr, ymirr, zmirr;
  double rcoordinate, rotation;
  MapCoordinates(x, y, z, xmirr, ymirr, zmirr, rcoordinate, rotation);

  // Do not proceed if not properly initialised.
  if (!m_ready) {
    PrintNotReady("GetMedium");
    return nullptr;
  }
  if (m_warning) PrintWarning("GetMedium");

  // Find the element that contains this point
  double t1, t2, t3, t4, jac[4][4], det;
  const int imap = FindElement13(x, y, z, t1, t2, t3, t4, jac, det);
  if (imap < 0) {
    if (m_debug) {
      std::cout << m_className << "::GetMedium:\n";
      std::cout << "    Point (" << x << ", " << y << ", " << z
                << ") not in the mesh.\n";
    }
    return nullptr;
  }
  const Element& element = elements[imap];
  if (element.matmap >= m_nMaterials) {
    if (m_debug) {
      std::cerr << m_className << "::GetMedium:\n";
      std::cerr << "    Point (" << x << ", " << y
                << ") has out of range material number " << imap << ".\n";
    }
    return nullptr;
  }

  if (m_debug) {
    PrintElement("GetMedium", x, y, z, t1, t2, t3, t4, element, 10);
  }

  return materials[element.matmap].medium;
}

double ComponentComsol::GetElementVolume(const unsigned int i) {
  if (i >= elements.size()) return 0.;
  const Element& element = elements[i];
  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];

  // Uses formula V = |a (dot) b x c|/6
  // with a => "3", b => "1", c => "2" and origin "0"
  const double vol =
      fabs((n3.x - n0.x) *
               ((n1.y - n0.y) * (n2.z - n0.z) - (n2.y - n0.y) * (n1.z - n0.z)) +
           (n3.y - n0.y) *
               ((n1.z - n0.z) * (n2.x - n0.x) - (n2.z - n0.z) * (n1.x - n0.x)) +
           (n3.z - n0.z) * ((n1.x - n0.x) * (n2.y - n0.y) -
                            (n3.x - n0.x) * (n1.y - n0.y))) /
      6.;
  return vol;
}

void ComponentComsol::GetAspectRatio(const unsigned int i, double& dmin,
                                     double& dmax) {
  if (i >= elements.size()) {
    dmin = dmax = 0.;
    return;
  }

  const Element& element = elements[i];
  const int np = 4;
  // Loop over all pairs of vertices.
  for (int j = 0; j < np - 1; ++j) {
    const Node& nj = nodes[element.emap[j]];
    for (int k = j + 1; k < np; ++k) {
      const Node& nk = nodes[element.emap[k]];
      // Compute distance.
      const double dx = nj.x - nk.x;
      const double dy = nj.y - nk.y;
      const double dz = nj.z - nk.z;
      const double dist = sqrt(dx * dx + dy * dy + dz * dz);
      if (k == 1) {
        dmin = dmax = dist;
      } else {
        if (dist < dmin) dmin = dist;
        if (dist > dmax) dmax = dist;
      }
    }
  }
}

}  // namespace Garfield
