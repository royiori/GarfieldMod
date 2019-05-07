#ifndef G_DRIFTLINE_RKF_H
#define G_DRIFTLINE_RKF_H

#include <string>
#include <vector>

#include "GeometryBase.hh"
#include "Medium.hh"
#include "Sensor.hh"
#include "ViewDrift.hh"

namespace Garfield {

/// Calculation of drift lines based on macroscopic transport coefficients
/// using Runge-Kutta-Fehlberg integration.

class DriftLineRKF {
 public:
  DriftLineRKF();
  ~DriftLineRKF() {}

  void SetSensor(Sensor* s);

  void EnablePlotting(ViewDrift* view);
  void DisablePlotting();

  void SetIntegrationAccuracy(const double a);
  void SetMaximumStepSize(const double ms);

  void EnableStepSizeLimit() { m_useStepSizeLimit = true; }
  void DisableStepSizeLimit() { m_useStepSizeLimit = false; }

  void SetMaxSteps(const unsigned int m) { m_maxSteps = m; }
  void EnableRejectKinks() { m_rejectKinks = true; }
  void DisableRejectKinks() { m_rejectKinks = false; }

  void SetElectronSignalScalingFactor(const double scale) {
    m_scaleElectronSignal = scale;
  }
  void SetHoleSignalScalingFactor(const double scale) {
    m_scaleHoleSignal = scale;
  }
  void SetIonSignalScalingFactor(const double scale) {
    m_scaleIonSignal = scale;
  }

  bool DriftElectron(const double x0, const double y0, const double z0,
                     const double t0);
  bool DriftHole(const double x0, const double y0, const double z0,
                 const double t0);
  bool DriftIon(const double x0, const double y0, const double z0,
                const double t0);

  void GetEndPoint(double& x, double& y, double& z, double& t, int& st) const;
  unsigned int GetNumberOfDriftLinePoints() const { return m_nPoints; }
  void GetDriftLinePoint(const unsigned int i, double& x, double& y, double& z,
                         double& t) const;

  double GetArrivalTimeSpread();
  double GetGain();
  double GetDriftTime() const {
    return m_nPoints > 0 ? m_path[m_nPoints - 1].t : 0.;
  }

  void EnableDebugging() { m_debug = true; }
  void DisableDebugging() { m_debug = false; }

  void EnableVerbose() { m_verbose = true; }
  void DisableVerbose() { m_verbose = false; }

 private:
  static const unsigned int ParticleTypeElectron = 0;
  static const unsigned int ParticleTypeIon = 1;
  static const unsigned int ParticleTypeHole = 2;

  std::string m_className = "DriftLineRKF";

  Sensor* m_sensor = nullptr;
  Medium* m_medium = nullptr;

  unsigned int m_particleType;
  double m_maxStepSize = 1.e8;
  double m_accuracy = 1.e-8;
  unsigned int m_maxSteps = 1000;
  unsigned int m_maxStepsToWire = 1000;
  bool m_rejectKinks = true;
  bool m_useStepSizeLimit = false;

  ViewDrift* m_view = nullptr;

  struct step {
    // Position
    double x;
    double y;
    double z;
    // Time
    double t;
    // Integrated Townsend coefficient
    double alphaint;
  };
  std::vector<step> m_path;
  int m_status = 0;
  unsigned int m_nPoints = 0;

  double m_scaleElectronSignal = 1.;
  double m_scaleHoleSignal = 1.;
  double m_scaleIonSignal = 1.;

  bool m_debug = false;
  bool m_verbose = false;

  // Calculate a drift line starting at a given position.
  bool DriftLine(const double x0, const double y0, const double z0,
                 const double t0);
  // Calculate transport parameters for the respective particle type.
  bool GetVelocity(const double x, const double y, const double z, double& vx,
                   double& vy, double& vz, int& status);
  bool GetVelocity(const double ex, const double ey, const double ez,
                   const double bx, const double by, const double bz,
                   double& vx, double& vy, double& vz) const;
  bool GetDiffusion(const double ex, const double ey, const double ez,
                    const double bx, const double by, const double bz,
                    double& dl, double& dt) const;
  bool GetTownsend(const double ex, const double ey, const double ez,
                   const double bx, const double by, const double bz,
                   double& alpha) const;
  // Terminate a drift line at the edge of a boundary.
  bool EndDriftLine();
  // Drift a particle to a wire
  bool DriftToWire(const double xw, const double yw, const double rw);
  // Determine the longitudinal diffusion over the drift line.
  double IntegrateDiffusion(const double x, const double y, const double z,
                            const double xe, const double ye, const double ze);
  // Determine the effective gain over the drift line.
  double IntegrateTownsend(const double x, const double y, const double z,
                           const double xe, const double ye, const double ze,
                           const double tol);
  // Calculate the signal for the current drift line.
  void ComputeSignal(const double q, const double scale) const;
};
}

#endif
