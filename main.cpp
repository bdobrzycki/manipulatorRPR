/*
* Bartosz Dobrzycki, 4 November 2013
*
* The following example demonstrates how to solve
* Revolute-Prismatic-Revolute (RPR) robotic manipulator forward kinematics
* problem using quaternions, please see the attached image: manipulatorRPR.pdf
*
* NOTE: Assertions and some parts of the functionality have been removed
* from the code for clarity (e.g -= operator or Magnitude() for Vector3 etc.).
* The code represents the bare minimum implementation-wise to solve the problem.
*
* Please see the attached 'QuaternionRotationOperator.pdf' paper that demonstrates
* how to derive an efficient algorithm (used here) that rotates input vector with a quaternion,
* utilizing the special rule for quaternion multiplication.
*/
#include <gl\glut.h>
#include <cmath>
#include <vector>
#include <iostream>

namespace MathTools
{
  // Define constant values.
  const double PI = 3.14159265358979323846264338;
  const double PI_OVER_TWO = 1.57079632679489661923132169;
  const double PI_OVER_FOUR = 0.78539816339744830961566084;
  const double ONE_EIGHTY_OVER_PI = 57.29577951308232087679815481;

  // Function converts angle in radians to angle in degrees.
  inline double Rad2Deg(double rad) { return rad * ONE_EIGHTY_OVER_PI; }

  //////////////////////////////////////////////////////////////////////
  // Function template that returns minimum of two values.
  template <typename T>
  inline const T& minimum(const T& a, const T& b)
  {
    return (a < b) ? a : b;
  }
  ////////////////////////////////////////////////////////////////////
  // Function template that returns maximum of two values.
  template <typename T>
  inline const T& maximum(const T& a, const T& b)
  {
    return (a < b) ? b : a;
  }
}
//////////////////////////////////////////////////////////////////////////
class Quaternion;
//////////////////////////////////////////////////////////////////////
// Three dimensional vector class.
class Vector3
{
friend
  class Quaternion; // Allows Quaternion class to access the private
                    // members of Vector3 class.
private:
  double x, y, z;   // Holds vector x, y, z components.

public:
  // Explicit default ctor.
  explicit Vector3(double _x = 0.0, double _y = 0.0, double _z = 0.0)
    :x(_x), y(_y), z(_z)
  {}

  // Getter functions for vector components.
  inline double GetX() const { return x; }
  inline double GetY() const { return y; }
  inline double GetZ() const { return z; }

  // Setter function that sets vector x, y, z components.
  void SetXYZ(double _x, double _y, double _z)
  {
    x = _x; y = _y; z = _z;
  }

  inline const Vector3 operator+(const Vector3& w) const
  {
    return Vector3(x + w.x, y + w.y, z + w.z);
  }

  Vector3& operator+=(const Vector3& w)
  {
    x += w.x; y += w.y; z += w.z;
    return *this;
  }

  inline const Vector3 operator-(const Vector3& w) const
  {
    return Vector3(x - w.x, y - w.y, z - w.z);
  }

  // Vector multiplication by a scalar operator.
  inline const Vector3 operator*(double a) const
  {
    return Vector3(x * a, y * a, z * a);
  }

  // Vector scalar product operator.
  inline const double operator*(const Vector3& v) const
  {
    return x * v.x + y * v.y + z * v.z;
  }

  // Vector cross product operator.
  inline const Vector3 operator%(const Vector3& v) const
  {
    return Vector3(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
  }
};
//////////////////////////////////////////////////////////////////////
// In this implementation Quaternion is represented by a real part r and
// three dimensional vector part v.
class Quaternion
{
private:
  double r;  //< Holds real part of the quaternion.
  Vector3 v; //< Holds vector part of the quaternion.

public:
  // Default ctor creates identity quaternion.
  Quaternion()
    : r(1.0), v(0.0, 0.0, 0.0)
  {}

  // Explicit ctor creates quaternion from real and vector parts.
  explicit Quaternion(double _r, const Vector3& _v)
    : r(_r), v(_v)
  {}

  // Function returns the real part of the quaternion.
  inline double GetRealPart() const { return r; }

  // Function returns the vector part of the quaternion.
  inline Vector3 GetVectorPart() const { return v; }

  // Function creates quaternion form unit axis and angle (in radians).
  inline void FromAxisAngle(const Vector3& axis, double angle)
  {
    angle /= 2.0;
    r = cos(angle);
    v = axis * sin(angle);
  }

  // Function returns angle in radians.
  inline double GetAngle() const
  {
    // Make sure acos() function parameter is in [-1.0, 1.0] range.
    const double realClamped =
      (r > 0.0) ? MathTools::minimum(r, 1.0) : MathTools::maximum(r, -1.0);
    return 2.0 * acos(realClamped);
  }

  // Set quaternion with new real and vector parts.
  inline void Set(double _r, const Vector3& _v)
  {
    r = _r;
    v = _v;
  }

  // Multiplication operator which follows exact mathematical definition but
  // becomes inefficient due to excessive usage of dot and cross product vector operators.
  // s = qp = qrpr - q · p + qrp + prq + q × p
  inline const Quaternion operator^(const Quaternion& p) const
  {
    return Quaternion(
      r*p.r - v*p.v,                   //< qrpr - q · p
      Vector3(p.v*r + v*p.r + v%p.v)); //< qrp + prq + q × p
  }

  // Multiplication operator in its unwrapped and efficient implementation.
  const Quaternion operator*(const Quaternion& p) const
  {
    return Quaternion(
      r*p.r - v.x*p.v.x-v.y*p.v.y - v.z*p.v.z,
      Vector3(r*p.v.x + p.r*v.x + v.y*p.v.z - v.z*p.v.y,
              r*p.v.y - v.x*p.v.z + p.r*v.y + v.z*p.v.x,
              r*p.v.z + v.x*p.v.y - v.y*p.v.x + p.r*v.z));
  }

  // Complex conjugate postfix operator.
  inline const Quaternion operator~() const
  {
    return Quaternion(r, v*(-1.0));
  }

  // Function makes this quaternion a unit length quaternion.
  void Normalise()
  {
    const double length = sqrt(r*r + v.x*v.x + v.y*v.y + v.z*v.z);
    if (length)
    {
      r /= length;
      v = v * (1.0 / length);
    }
  }

  // Function rotates input vector w with a quaternion i.e. it explicitly applies
  // quaternion rotation operator to that vector.
  // It returns new, rotated vector, extracted from the pure quaternion.
  // Implementation strictly follows mathematical definition, but it is computationally inefficient.
  inline const Vector3 Rotate(const Vector3& w) const
  {
    // Assuming e has unit length here.
    const Quaternion& e = *this;
    return (e*(Quaternion(0.0, w)*(~e))).v;
  }

  // Function rotates input vector with a quaternion in most efficient way.
  // Please see attached 'QuaternionRotationOperator.pdf' paper to see
  // how to derive the efficient formula.
  // w' = e(0+w)e* = (2*r*r - 1)w + 2r(v × w) + 2(v · w)v
  inline const Vector3 RotateFast(const Vector3& w) const
  {
    // Assuming e has unit length here.
    const double exrMult = 2.0 * r;
    const double rMult = exrMult * r - 1.0;
    const double eMult = 2.0 * (v.x*w.x + v.y*w.y + v.z*w.z);
    return Vector3(
      rMult*w.x + exrMult*(v.y*w.z - v.z*w.y) + eMult*v.x,
      rMult*w.y + exrMult*(v.z*w.x - v.x*w.z) + eMult*v.y,
      rMult*w.z + exrMult*(v.x*w.y - v.y*w.x) + eMult*v.z);
  }
};
//////////////////////////////////////////////////////////////////////
// Color line for debug drawing.
class GLLine
{
private:
  GLfloat from[3];
  GLfloat to[3];
  GLfloat color[3];

public:
  explicit GLLine(
    const Vector3& from,
    const Vector3& to,
    const Vector3& color)
  {
    GLLine* const line = this;
    line->from[0] = GLfloat(from.GetX());
    line->from[1] = GLfloat(from.GetY());
    line->from[2] = GLfloat(from.GetZ());

    line->to[0] = GLfloat(to.GetX());
    line->to[1] = GLfloat(to.GetY());
    line->to[2] = GLfloat(to.GetZ());

    line->color[0] = GLfloat(color.GetX());
    line->color[1] = GLfloat(color.GetY());
    line->color[2] = GLfloat(color.GetZ());
  }

  void Draw() const
  {
    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_LINES);
    glVertex3f(from[0], from[1], from[2]);
    glVertex3f(to[0], to[1], to[2]);
    glEnd();
  }
};
//////////////////////////////////////////////////////////////////////
// Define constant values.
const double timeDelta = 0.01;
const Vector3 zero(0.0, 0.0, 0.0);

// Global reference frame G axes.
const Vector3 xAxis (1.0, 0.0, 0.0);
const Vector3 yAxis (0.0, 1.0, 0.0);
const Vector3 zAxis (0.0, 0.0, 1.0);

// Colors.
const Vector3 red(1.0, 0.0, 0.0);
const Vector3 green(0.0, 1.0, 0.0);
const Vector3 blue(0.0, 0.0, 1.0);
const Vector3 orange(1.0, 0.5, 0.0);
const Vector3 white(1.0, 1.0, 1.0);

// Buffer for storing lines for drawing.
std::vector<GLLine> lines;
// Line buffer constant iterator.
typedef std::vector<GLLine>::const_iterator GLLineConstIter;

//////////////////////////////////////////////////////////////////////
// Simulation step function.
void SimStep (double timeDelta)
{
  static double time = 0.0;
  time += timeDelta;
 
  // Manipulator setup;
  // please see the attached manipulatorRPR.pdf diagram.

  // B2rP is the position of the tip of the manipulator P
  // expressed in local reference frame B2.
  const Vector3 B2rP(0.0, 2.0, 0.0);

  // Reference frame B2 can rotate about its local x2 axis,
  // relatively to B1 frame.
  // This models a single degree of freedom revolute joint,
  // that can rotate in range: <Pi/4, -Pi/4>.
  const double phi2 = (MathTools::PI_OVER_FOUR) * sin(time);  //< Rotation angle.
  const Vector3 u2(1.0, 0.0, 0.0);  //< Rotation axis (local x2).
  // Quaternion representing a rotation of B2 frame relatively to B1 frame.
  Quaternion B1eB2;
  B1eB2.FromAxisAngle(u2, phi2);

  // Prismatic joint.
  // Position of B2 frame expressed in local B1 frame.
  const Vector3 B1dB2(0.0, 4.0*abs(sin(time)), 0.0);

  // Reference frame B1 can rotate about global Z axis,
  // relatively to the global reference frame G.
  // This models a single degree of freedom revolute joint,
  // that can rotate in range: <Pi, -Pi>.
  const double phi1 = MathTools::PI * sin(time);
  const Vector3 u1(0.0, 0.0, 1.0);  //< Rotation axis (Global Z).
  // Quaternion representing a rotation of B1 frame
  // relatively to the global reference frame G.
  Quaternion GeB1;
  GeB1.FromAxisAngle(u1, phi1);

  // Position of frame B1 relatively to the global reference frame G.
  const Vector3 GdB1(0.0, 0.0, 2.0);

  // Calculate rotation of B2 frame relatively to global reference frame G.
  const Quaternion GeB2 = GeB1 * B1eB2;

  ////////////////////////////////////////////////////////////////////
  // Solve the FK utilizing the quaternion special product only.
  // This is strict mathematically, but inefficient.
  /*const Vector3 GrP =
    (GeB2*(Quaternion(0.0, B2rP)*(~GeB2))).GetVectorPart() +
    (GeB1*(Quaternion(0.0, B1dB2)*(~GeB1))).GetVectorPart() +
    GdB1;*/

  ////////////////////////////////////////////////////////////////////
  // Solve the FK utilizing the quaternion special product only,
  // this time use more human friendly multiplication formula, defined
  // as ^ operator.
  // This is strict mathematically, but inefficient.
  /*const Vector3 GrP =
    (GeB2^(Quaternion(0.0, B2rP)^(~GeB2))).GetVectorPart() +
    (GeB1^(Quaternion(0.0, B1dB2)^(~GeB1))).GetVectorPart() +
    GdB1;*/

  ////////////////////////////////////////////////////////////////////
  // Solve the FK in the most efficient way.
  const Vector3 GrP = GeB2.RotateFast(B2rP) + GeB1.RotateFast(B1dB2) + GdB1;

  ////////////////////////////////////////////////////////////////////
  // Manipulator debug-draw in global reference frame G.
  // B1 reference frame.
  lines.push_back(GLLine(GdB1, GdB1 + GeB1.RotateFast(xAxis), red));
  lines.push_back(GLLine(GdB1, GdB1 + GeB1.RotateFast(yAxis), green));
  lines.push_back(GLLine(GdB1, GdB1 + GeB1.RotateFast(zAxis), blue));

  // B2 reference frame translation relatively to B1 reference frame.
  const Vector3 GB2 = GdB1 + GeB1.RotateFast(B1dB2);
  lines.push_back(GLLine(GdB1, GB2, orange));

  // B2 reference frame.
  lines.push_back(GLLine(GB2, GB2 + GeB2.RotateFast(xAxis), red));
  lines.push_back(GLLine(GB2, GB2 + GeB2.RotateFast(yAxis), green));
  lines.push_back(GLLine(GB2, GB2 + GeB2.RotateFast(zAxis), blue));

  // Tip of the manipulator P.
  lines.push_back(GLLine(GB2, GrP, white));

  // Print out joint current angles.
  std::cout << "RevJ1: angle(DEG): " << MathTools::Rad2Deg(GeB1.GetAngle())
            << " RevJ2: angle(DEG): " << MathTools::Rad2Deg(B1eB2.GetAngle()) << std::endl;
}
//////////////////////////////////////////////////////////////////////
void Init(void)
{
  glShadeModel(GL_SMOOTH);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f, (GLfloat)400/(GLfloat)400, 1.0f, 1000.0f);
  glEnable(GL_DEPTH_TEST);
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
}
//////////////////////////////////////////////////////////////////////
void Display(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(10.0, 10.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

  lines.push_back(GLLine(zero, xAxis, Vector3(1, 0, 0)));
  lines.push_back(GLLine(zero, yAxis, Vector3(0, 1, 0)));
  lines.push_back(GLLine(zero, zAxis, Vector3(0, 0, 1)));

  // Draw lines.
  GLLineConstIter iter = lines.begin();
  for (iter; iter != lines.end(); ++iter)
  {
    iter->Draw();
  }

  lines.clear();

  glutSwapBuffers();
}
//////////////////////////////////////////////////////////////////////
void TimerFunction(int value)
{
  SimStep(timeDelta);
  Display();
  glutTimerFunc(10, TimerFunction, 1);
}
//////////////////////////////////////////////////////////////////////
int main (int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(400, 400);
  glutInitWindowPosition(380, 180);
  glutCreateWindow("Revolute-Prismatic-Revolute (RPR) Manipulator");
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glutDisplayFunc(Display);
  glutTimerFunc(10, TimerFunction, 1);
  Init();
  glutMainLoop();
  std::cin.get();
  return 0;
}
