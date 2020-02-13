package LMC_MBS
  function waveNumberIterator
    input Real d "Water depth";
    input Real omega[:] "Wave frequency";
    output Real k[size(omega, 1)] "Wave number";
  protected
    constant Real g = Modelica.Constants.g_n;
    constant Real pi = Modelica.Constants.pi;
    parameter Integer n = size(omega, 1);
    Real T[size(omega, 1)] "Wave period";
    Real L0[size(omega, 1)] "Deepwater wavelength";
    Real Ll(start = 0, fixed = true) "Intermediate loop value";
    Real Llc(start = 0, fixed = true) "Intermediate loop comparator value";
    Real L[size(omega, 1)] "Iterated wave length";
  algorithm
    T := 2 * pi ./ omega;
    L0 := g * T .^ 2 / (2 * pi);
    for i in 1:size(omega, 1) loop
      Ll := L0[i];
      Llc := 0;
      while abs(Llc - Ll) > 0.001 loop
        Llc := Ll;
        L[i] := g * T[i] ^ 2 / (2 * pi) * tanh(2 * pi / Ll * d);
        Ll := L[i];
      end while;
    end for;
    k := 2 * pi ./ L;
    annotation(
      Diagram(coordinateSystem(initialScale = 0.05)));
  end waveNumberIterator;

  function linearInterpolatorSV
    input Real x[:];
    input Real y[:];
    input Real p;
    output Real q;
  protected
    Integer j;
  algorithm
    if p == x[size(x, 1)] then
      q := y[size(x, 1)];
    else
      j := 1;
      while p >= x[j + 1] loop
        j := j + 1;
      end while;
      q := (y[j + 1] - y[j]) / (x[j + 1] - x[j]) * (p - x[j]) + y[j];
    end if;
  end linearInterpolatorSV;

  model UP_Seg
    import Modelica.Mechanics.MultiBody.Types;
    import Modelica.SIunits.Conversions.to_unit1;
    import SI = Modelica.SIunits;
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a "Coordinate system fixed to the component with one cut-force and cut-torque" annotation(
      Placement(transformation(extent = {{-116, -16}, {-84, 16}})));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b "Coordinate system fixed to the component with one cut-force and cut-torque" annotation(
      Placement(transformation(extent = {{84, -16}, {116, 16}})));
    Modelica.Blocks.Interfaces.RealOutput y[3] "Connector of Real output signals" annotation(
      Placement(transformation(extent = {{-80, -80}, {-60, -60}}, rotation = 270)));
    parameter Boolean animation = true "= true, if animation shall be enabled";
    parameter SI.Position r[3](start = {0, 0, 0}) "Vector from frame_a to frame_b resolved in frame_a";
    parameter Types.ShapeType shapeType = "cylinder" "Type of shape" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter SI.Position r_shape[3] = {0, 0, 0} "Vector from frame_a to shape origin, resolved in frame_a" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter Types.Axis lengthDirection = to_unit1(r - r_shape) "Vector in length direction of shape, resolved in frame_a" annotation(
      Evaluate = true,
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter Types.Axis widthDirection = {0, 1, 0} "Vector in width direction of shape, resolved in frame_a" annotation(
      Evaluate = true,
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter SI.Length length = Modelica.Math.Vectors.length(r - r_shape) "Length of shape" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter SI.Distance width = length / world.defaultWidthFraction "Width of shape" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter SI.Distance height = width "Height of shape" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter Types.ShapeExtra extra = 0.0 "Additional parameter depending on shapeType (see docu of Visualizers.Advanced.Shape)" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    input Types.Color color = Modelica.Mechanics.MultiBody.Types.Defaults.RodColor "Color of shape" annotation(
      Dialog(colorSelector = true, tab = "Animation", group = "if animation = true", enable = animation));
    input Types.SpecularCoefficient specularCoefficient = world.defaultSpecularCoefficient "Reflection of ambient light (= 0: light is completely absorbed)" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
  protected
    outer Modelica.Mechanics.MultiBody.World world;
    Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape shape(shapeType = shapeType, color = color, specularCoefficient = specularCoefficient, r_shape = r_shape, lengthDirection = lengthDirection, widthDirection = widthDirection, length = length, width = width, height = height, extra = extra, r = frame_a.r_0, R = frame_a.R) if world.enableAnimation and animation;
  equation
    Connections.branch(frame_a.R, frame_b.R);
    assert(cardinality(frame_a) > 0 or cardinality(frame_b) > 0, "Neither connector frame_a nor frame_b of FixedTranslation object is connected");
    frame_b.r_0 = frame_a.r_0 + Modelica.Mechanics.MultiBody.Frames.resolve1(frame_a.R, r);
    frame_b.R = frame_a.R;
/* Force and torque balance */
    zeros(3) = frame_a.f + frame_b.f;
    zeros(3) = frame_a.t + frame_b.t + cross(r, frame_b.f);
    y[1] = frame_a.r_0[1];
    y[2] = frame_a.r_0[2];
    y[3] = frame_a.r_0[3];
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-99, 5}, {101, -5}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 85}, {150, 45}}, textString = "%name", lineColor = {0, 0, 255}), Text(extent = {{150, -50}, {-150, -20}}, lineColor = {0, 0, 0}, textString = "r=%r"), Text(extent = {{-89, 38}, {-53, 13}}, lineColor = {128, 128, 128}, textString = "a"), Text(extent = {{57, 39}, {93, 14}}, lineColor = {128, 128, 128}, textString = "b")}),
      Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-100, 5}, {100, -5}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-95, 20}, {-58, 20}}, color = {128, 128, 128}, arrow = {Arrow.None, Arrow.Filled}), Line(points = {{-94, 18}, {-94, 50}}, color = {128, 128, 128}, arrow = {Arrow.None, Arrow.Filled}), Text(extent = {{-72, 35}, {-58, 24}}, lineColor = {128, 128, 128}, textString = "x"), Text(extent = {{-113, 57}, {-98, 45}}, lineColor = {128, 128, 128}, textString = "y"), Line(points = {{-100, -4}, {-100, -69}}, color = {128, 128, 128}), Line(points = {{-100, -63}, {90, -63}}, color = {128, 128, 128}), Text(extent = {{-22, -39}, {16, -63}}, lineColor = {128, 128, 128}, textString = "r"), Polygon(points = {{88, -59}, {88, -68}, {100, -63}, {88, -59}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{100, -3}, {100, -68}}, color = {128, 128, 128}), Line(points = {{69, 20}, {106, 20}}, color = {128, 128, 128}, arrow = {Arrow.None, Arrow.Filled}), Line(points = {{70, 18}, {70, 50}}, color = {128, 128, 128}, arrow = {Arrow.None, Arrow.Filled}), Text(extent = {{92, 35}, {106, 24}}, lineColor = {128, 128, 128}, textString = "x"), Text(extent = {{51, 57}, {66, 45}}, lineColor = {128, 128, 128}, textString = "y")}),
      Documentation(info = "<html>
  <p>
  To be added.
  <html>"));
  end UP_Seg;

  model LO_Seg
    import Modelica.Mechanics.MultiBody.Types;
    import Modelica.SIunits.Conversions.to_unit1;
    import SI = Modelica.SIunits;
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a "Coordinate system fixed to the component with one cut-force and cut-torque" annotation(
      Placement(transformation(extent = {{-116, -16}, {-84, 16}})));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b "Coordinate system fixed to the component with one cut-force and cut-torque" annotation(
      Placement(transformation(extent = {{84, -16}, {116, 16}})));
    Modelica.Blocks.Interfaces.RealOutput y[3] "Connector of Real output signals" annotation(
      Placement(transformation(extent = {{60, -80}, {80, -60}}, rotation = 270)));
    parameter Boolean animation = true "= true, if animation shall be enabled";
    parameter SI.Position r[3](start = {0, 0, 0}) "Vector from frame_a to frame_b resolved in frame_a";
    parameter Types.ShapeType shapeType = "cylinder" "Type of shape" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter SI.Position r_shape[3] = {0, 0, 0} "Vector from frame_a to shape origin, resolved in frame_a" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter Types.Axis lengthDirection = to_unit1(r - r_shape) "Vector in length direction of shape, resolved in frame_a" annotation(
      Evaluate = true,
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter Types.Axis widthDirection = {0, 1, 0} "Vector in width direction of shape, resolved in frame_a" annotation(
      Evaluate = true,
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter SI.Length length = Modelica.Math.Vectors.length(r - r_shape) "Length of shape" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter SI.Distance width = length / world.defaultWidthFraction "Width of shape" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter SI.Distance height = width "Height of shape" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    parameter Types.ShapeExtra extra = 0.0 "Additional parameter depending on shapeType (see docu of Visualizers.Advanced.Shape)" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
    input Types.Color color = Modelica.Mechanics.MultiBody.Types.Defaults.RodColor "Color of shape" annotation(
      Dialog(colorSelector = true, tab = "Animation", group = "if animation = true", enable = animation));
    input Types.SpecularCoefficient specularCoefficient = world.defaultSpecularCoefficient "Reflection of ambient light (= 0: light is completely absorbed)" annotation(
      Dialog(tab = "Animation", group = "if animation = true", enable = animation));
  protected
    outer Modelica.Mechanics.MultiBody.World world;
    Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape shape(shapeType = shapeType, color = color, specularCoefficient = specularCoefficient, r_shape = r_shape, lengthDirection = lengthDirection, widthDirection = widthDirection, length = length, width = width, height = height, extra = extra, r = frame_a.r_0, R = frame_a.R) if world.enableAnimation and animation;
  equation
    Connections.branch(frame_a.R, frame_b.R);
    assert(cardinality(frame_a) > 0 or cardinality(frame_b) > 0, "Neither connector frame_a nor frame_b of FixedTranslation object is connected");
    frame_b.r_0 = frame_a.r_0 + Modelica.Mechanics.MultiBody.Frames.resolve1(frame_a.R, r);
    frame_b.R = frame_a.R;
/* Force and torque balance */
    zeros(3) = frame_a.f + frame_b.f;
    zeros(3) = frame_a.t + frame_b.t + cross(r, frame_b.f);
    y[1] = frame_b.r_0[1];
    y[2] = frame_b.r_0[2];
    y[3] = frame_b.r_0[3];
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-99, 5}, {101, -5}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 85}, {150, 45}}, textString = "%name", lineColor = {0, 0, 255}), Text(extent = {{150, -50}, {-150, -20}}, lineColor = {0, 0, 0}, textString = "r=%r"), Text(extent = {{-89, 38}, {-53, 13}}, lineColor = {128, 128, 128}, textString = "a"), Text(extent = {{57, 39}, {93, 14}}, lineColor = {128, 128, 128}, textString = "b")}),
      Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-100, 5}, {100, -5}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{-95, 20}, {-58, 20}}, color = {128, 128, 128}, arrow = {Arrow.None, Arrow.Filled}), Line(points = {{-94, 18}, {-94, 50}}, color = {128, 128, 128}, arrow = {Arrow.None, Arrow.Filled}), Text(extent = {{-72, 35}, {-58, 24}}, lineColor = {128, 128, 128}, textString = "x"), Text(extent = {{-113, 57}, {-98, 45}}, lineColor = {128, 128, 128}, textString = "y"), Line(points = {{-100, -4}, {-100, -69}}, color = {128, 128, 128}), Line(points = {{-100, -63}, {90, -63}}, color = {128, 128, 128}), Text(extent = {{-22, -39}, {16, -63}}, lineColor = {128, 128, 128}, textString = "r"), Polygon(points = {{88, -59}, {88, -68}, {100, -63}, {88, -59}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Line(points = {{100, -3}, {100, -68}}, color = {128, 128, 128}), Line(points = {{69, 20}, {106, 20}}, color = {128, 128, 128}, arrow = {Arrow.None, Arrow.Filled}), Line(points = {{70, 18}, {70, 50}}, color = {128, 128, 128}, arrow = {Arrow.None, Arrow.Filled}), Text(extent = {{92, 35}, {106, 24}}, lineColor = {128, 128, 128}, textString = "x"), Text(extent = {{51, 57}, {66, 45}}, lineColor = {128, 128, 128}, textString = "y")}),
      Documentation(info = "<html>
  <p>
  To be added.
  <html>"));
  end LO_Seg;

  block DnB
    extends Modelica.Blocks.Icons.Block;
    Modelica.Blocks.Interfaces.RealInput c2[3] "Connector of Real input signals" annotation(
      Placement(transformation(extent = {{-120, -10}, {-100, 10}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealInput c1[3] "Connector of Real input signals" annotation(
      Placement(transformation(extent = {{100, -10}, {120, 10}}, rotation = 180)));
    Modelica.Blocks.Interfaces.RealOutput y1[3] "Connector of Real output signal" annotation(
      Placement(transformation(extent = {{-10, -120}, {10, -100}}, rotation = -90)));
    // import constants
    parameter Real pi = Modelica.Constants.pi;
    parameter Real g = Modelica.Constants.g_n;
    // Specify Environment Data
    parameter Real rho_w = 1025 "Density of sea water [kg/m^3]";
    parameter Real d = 50 "Constant water depth [m]";
    parameter Real Hr = 2 "Height of regular wave [m]";
    parameter Real Tr = 10 "Period of regular wave [s]";
    parameter Real zcg[:] = {-50, -25, -10, 0} "z coordinates of the current profile";
    parameter Real Uf[:] = {0,0,1,2} "Fully developed current velocities at z_cg";
    // Specify Line Data
    parameter Real D = 0.04 "Hydraulic diameter";
    parameter Real Db = 0.04 "Buoyancy diameter";
    parameter Real Ca_n = 1 "Normal added mass coefficient";
    parameter Real Ca_t = 0.5 "Tangential added mass coefficient";
    parameter Real Cd_n = 1 "Normal drag coefficient";
    parameter Real Cd_t = 0.25 "Tangential drag coefficient";
    // Specify Ramping Data
    parameter Real Trmp = 10 "Ramp time for current and wave height [s]";
    // Calculate Environment Derived Parameters
    parameter Real omega[:] = {2 * pi / Tr} "Wave frequency [rad/s]";
    parameter Real k[:] = waveNumberIterator(d, omega) "Function call for wave number k [m^-1]";
    // Calculate Line Derived Parameters
    parameter Real Cm_n = 1 + Ca_n "Normal inertia coefficient";
    parameter Real Cm_t = 1 + Ca_t "Tangential inertia coefficient";
    // Declare simulation variables
    Real seg_l;
    Real stheta;
    Real ctheta;
    Real x_m;
    Real y_m;
    Real Ucg[size(Uf, 1)];
    Real SSE;
    Real c_m;
    Real c_n;
    Real c_t;
    Real u;
    Real w;
    Real v_wn;
    Real v_wt;
    Real du;
    Real dw;
    Real a_wn;
    Real a_wt;
    Real v_mx;
    Real v_my;
    Real v_mn;
    Real v_mt;
    Real dv_mx;
    Real dv_my;
    Real a_mn;
    Real a_mt;
    Real Mf_n;
    Real Mf_t;
    Real Mf_x;
    Real Mf_y;
    Real Bu;
  equation
// Calculate segment lengths and inclinations
    seg_l = sqrt((c1[1] - c2[1]) * (c1[1] - c2[1]) + (c1[2] - c2[2]) * (c1[2] - c2[2]) + (c1[3] - c2[3]) * (c1[3] - c2[3]));
    stheta = (c2[2] - c1[2]) / seg_l;
    ctheta = (c2[1] - c1[1]) / seg_l;
//Calculate segment mid-point instantaneous co-ordinates
    x_m = (c1[1] + c2[1]) / 2;
    y_m = (c1[2] + c2[2]) / 2;
//Determine instantaneous current profile (for ramping)
    for i in 1:size(Uf, 1) loop
      if time < Trmp then
        Ucg[i] = sin(pi/2*time / Trmp) * Uf[i];
      else
        Ucg[i] = Uf[i];
      end if;
    end for;
// Calculate sea-surface elevation directly above the segment mid-point
    if time < Trmp then
      SSE = sin(pi/2*time / Trmp) * Hr / 2 * cos(k[1] * x_m - omega[1] * time);
    else
      SSE = Hr / 2 * cos(k[1] * x_m - omega[1] * time);
    end if;
//Determine current velocity components at segment mid-point (current profile moved with SSE)
    c_m = linearInterpolatorSV(zcg, Ucg, y_m - SSE);
    c_n = c_m * stheta;
    c_t = c_m * ctheta;
//Determine wave-induced water-particle velocity and acceleration components at segment mid-point
    if time < Trmp then
      u = sin(pi/2*time / Trmp) * pi * Hr / Tr * cosh(k[1] * (y_m - SSE + d)) / sinh(k[1] * d) * cos(k[1] * x_m - omega[1] * time);
      w = sin(pi/2*time / Trmp) * pi * Hr / Tr * sinh(k[1] * (y_m - SSE + d)) / sinh(k[1] * d) * sin(k[1] * x_m - omega[1] * time);
      du = sin(pi/2*time / Trmp) * 2 * pi * pi * Hr / (Tr * Tr) * cosh(k[1] * (y_m - SSE + d)) / sinh(k[1] * d) * sin(k[1] * x_m - omega[1] * time);
      dw = -sin(pi/2*time / Trmp) * 2 * pi * pi * Hr / (Tr * Tr) * sinh(k[1] * (y_m - SSE + d)) / sinh(k[1] * d) * cos(k[1] * x_m - omega[1] * time);
    else
      u = pi * Hr / Tr * cosh(k[1] * (y_m - SSE + d)) / sinh(k[1] * d) * cos(k[1] * x_m - omega[1] * time);
      w = pi * Hr / Tr * sinh(k[1] * (y_m - SSE + d)) / sinh(k[1] * d) * sin(k[1] * x_m - omega[1] * time);
      du = 2 * pi * pi * Hr / (Tr * Tr) * cosh(k[1] * (y_m - SSE + d)) / sinh(k[1] * d) * sin(k[1] * x_m - omega[1] * time);
      dw = -2 * pi * pi * Hr / (Tr * Tr) * sinh(k[1] * (y_m - SSE + d)) / sinh(k[1] * d) * cos(k[1] * x_m - omega[1] * time);
    end if;
    v_wn = u * stheta - w * ctheta;
    v_wt = u * ctheta + w * stheta;
//du = der(u);
//dw = der(w);
    a_wn = du * stheta - dw * ctheta;
    a_wt = du * ctheta + dw * stheta;
// Determine segment mid point velocity and acceleration components
    v_mx = der(x_m);
    v_my = der(y_m);
    v_mn = v_mx * stheta - v_my * ctheta;
    v_mt = v_mx * ctheta + v_my * stheta;
    dv_mx = der(v_mx);
    dv_my = der(v_my);
    a_mn = dv_mx * stheta - dv_my * ctheta;
    a_mt = dv_mx * ctheta + dv_my * stheta;
//Determine drag components
    Mf_n = (Cm_n * rho_w * pi / 4 * D * D * a_wn - Ca_n * rho_w * pi / 4 * D * D * a_mn + 0.5 * Cd_n * rho_w * D * abs(v_wn + c_n - v_mn) * (v_wn + c_n - v_mn)) * seg_l;
    Mf_t = (Cm_t * rho_w * pi / 4 * D * D * a_wt - Ca_t * rho_w * pi / 4 * D * D * a_mt + 0.5 * Cd_t * rho_w * pi * D * abs(v_wt + c_t - v_mt) * (v_wt + c_t - v_mt)) * seg_l;
    Mf_x = Mf_n * stheta + Mf_t * ctheta;
    Mf_y = (-Mf_n * ctheta) + Mf_t * stheta;
    Bu = pi * Db * Db / 4 * seg_l * rho_w * g;
    y1[1] = Mf_x;
    y1[2] = Mf_y + Bu;
    y1[3] = 0;
    annotation(
      Documentation(info = "<html>
<p>
Block has a vector of continuous Real input signals and
one continuous Real output signal.
</p>
</html>"),
      experiment(StartTime = 0, StopTime = 20, Tolerance = 1e-06, Interval = 0.1),
      Icon(graphics = {Text(origin = {44, 32}, rotation = -90, extent = {{-42, 76}, {104, -154}}, textString = "D&B")}, coordinateSystem(initialScale = 0.1)));
  end DnB;
















  model SubseaLoad1
    Modelica.Mechanics.MultiBody.Parts.Fixed fixed1 annotation(
      Placement(visible = true, transformation(origin = {-70, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1 annotation(
      Placement(visible = true, transformation(origin = {-20, 80}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    LMC_MBS.UP_Seg uP_Seg1 annotation(
      Placement(visible = true, transformation(origin = {-20, 58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.PointMass LM1 annotation(
      Placement(visible = true, transformation(origin = {-20, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    LMC_MBS.LO_Seg lO_Seg1 annotation(
      Placement(visible = true, transformation(origin = {-20, 22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2 annotation(
      Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    LMC_MBS.UP_Seg uP_Seg2 annotation(
      Placement(visible = true, transformation(origin = {-20, -22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.PointMass LM2 annotation(
      Placement(visible = true, transformation(origin = {-20, -40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    LMC_MBS.LO_Seg lO_Seg2 annotation(
      Placement(visible = true, transformation(origin = {-20, -58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.PointMass Load annotation(
      Placement(visible = true, transformation(origin = {-20, -78}, extent = {{-16, -16}, {16, 16}}, rotation = -90)));
    LMC_MBS.DnB dnB1 annotation(
      Placement(visible = true, transformation(origin = {-70, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Forces.WorldForce force1 annotation(
      Placement(visible = true, transformation(origin = {-40, 40}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    LMC_MBS.DnB dnB2 annotation(
      Placement(visible = true, transformation(origin = {-70, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Forces.WorldForce force2 annotation(
      Placement(visible = true, transformation(origin = {-40, -40}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    inner Modelica.Mechanics.MultiBody.World world annotation(
      Placement(visible = true, transformation(origin = {-80, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(lO_Seg1.y, dnB1.c2) annotation(
      Line(points = {{-27, 15}, {-29, 15}, {-29, 19}, {-71, 19}, {-71, 27}}, color = {0, 0, 127}, thickness = 0.5));
    connect(force2.frame_b, LM2.frame_a) annotation(
      Line(points = {{-32, -40}, {-20, -40}, {-20, -40}, {-20, -40}}, color = {95, 95, 95}));
    connect(lO_Seg2.frame_b, Load.frame_a) annotation(
      Line(points = {{-20, -68}, {-20, -78}}, color = {95, 95, 95}));
    connect(LM2.frame_a, lO_Seg2.frame_a) annotation(
      Line(points = {{-20, -40}, {-20, -40}, {-20, -48}, {-20, -48}}, color = {95, 95, 95}));
    connect(lO_Seg2.y, dnB2.c2) annotation(
      Line(points = {{-25, -65}, {-27, -65}, {-27, -61}, {-69, -61}, {-69, -51}, {-69, -51}}, color = {0, 0, 127}, thickness = 0.5));
    connect(uP_Seg2.frame_b, LM2.frame_a) annotation(
      Line(points = {{-20, -32}, {-20, -32}, {-20, -40}, {-20, -40}}, color = {95, 95, 95}));
    connect(force1.frame_b, LM2.frame_a) annotation(
      Line(points = {{-32, 40}, {-20, 40}, {-20, 40}, {-20, 40}}, color = {95, 95, 95}));
    connect(uP_Seg2.y, dnB2.c1) annotation(
      Line(points = {{-27, -15}, {-27, -15}, {-27, -21}, {-71, -21}, {-71, -29}, {-71, -29}}, color = {0, 0, 127}, thickness = 0.5));
    connect(dnB2.y1, force1.force) annotation(
      Line(points = {{-57, -40}, {-49, -40}, {-49, -40}, {-49, -40}}, color = {0, 0, 127}, thickness = 0.5));
    connect(revolute2.frame_b, uP_Seg2.frame_a) annotation(
      Line(points = {{-12, 0}, {-6, 0}, {-6, -12}, {-20, -12}, {-20, -12}}, color = {95, 95, 95}));
    connect(revolute1.frame_b, uP_Seg1.frame_a) annotation(
      Line(points = {{-10, 80}, {-4, 80}, {-4, 68}, {-18, 68}, {-18, 68}}, color = {95, 95, 95}));
    connect(fixed1.frame_b, revolute1.frame_a) annotation(
      Line(points = {{-60, 80}, {-28, 80}}, color = {95, 95, 95}));
    connect(LM1.frame_a, lO_Seg1.frame_a) annotation(
      Line(points = {{-20, 40}, {-20, 40}, {-20, 32}, {-20, 32}}, color = {95, 95, 95}));
    connect(lO_Seg1.frame_b, revolute2.frame_a) annotation(
      Line(points = {{-18, 12}, {-32, 12}, {-32, 0}, {-26, 0}, {-26, 0}}, color = {95, 95, 95}));
    connect(uP_Seg1.frame_b, LM1.frame_a) annotation(
      Line(points = {{-18, 48}, {-18, 48}, {-18, 40}, {-18, 40}}, color = {95, 95, 95}));
    connect(uP_Seg1.y, dnB1.c1) annotation(
      Line(points = {{-25, 65}, {-26, 65}, {-26, 60}, {-68, 60}, {-68, 52}}, color = {0, 0, 127}, thickness = 0.5));
    connect(force1.frame_b, LM1.frame_a) annotation(
      Line(points = {{-32, 40}, {-20, 40}, {-20, 40}, {-20, 40}}, color = {95, 95, 95}));
    connect(dnB1.y1, force1.force) annotation(
      Line(points = {{-57, 40}, {-49, 40}, {-49, 40}, {-49, 40}}, color = {0, 0, 127}, thickness = 0.5));
    connect(force.frame_b, LM1.frame_a) annotation(
      Line(points = {{-12, 40}, {0, 40}, {0, 40}, {0, 40}}, color = {95, 95, 95}));
    connect(dnB1.y1, force.force) annotation(
      Line(points = {{-38, 40}, {-30, 40}, {-30, 40}, {-30, 40}}, color = {0, 0, 127}, thickness = 0.5));
  end SubseaLoad1;

  model SubseaLoad
    inner Modelica.Mechanics.MultiBody.World world annotation(
      Placement(visible = true, transformation(origin = {-70, -88}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.Fixed fixed1(r = {0, -2.5, 0}, r_shape = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {-70, 90}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(phi(fixed = true), useAxisFlange = true, w(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {-50, 90}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.UP_Seg uP_Seg1(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {-49, 73}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.PointMass LM1(m = 50, r_0(fixed = true, start = {0, -5.0, 0}), v_0(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {-49, 61}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    LMC_MBS.LO_Seg lO_Seg1(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {-49, 49}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    LMC_MBS.DnB dnB1 annotation(
      Placement(visible = true, transformation(origin = {-86, 60}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Forces.WorldForce force annotation(
      Placement(visible = true, transformation(origin = {-64, 60}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.DnB dnB2 annotation(
      Placement(visible = true, transformation(origin = {-86, 2}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Forces.WorldForce force1 annotation(
      Placement(visible = true, transformation(origin = {-64, 2}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.UP_Seg uP_Seg2(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {-49, 15}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.PointMass LM2(m = 50, r_0(fixed = true, start = {0, -10.0, 0}), v_0(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {-49, 3}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    LMC_MBS.LO_Seg lO_Seg2(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {-49, -9}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(phi(fixed = true), useAxisFlange = true, w(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {-50, 32}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.DnB dnB3 annotation(
      Placement(visible = true, transformation(origin = {-86, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Forces.WorldForce force2 annotation(
      Placement(visible = true, transformation(origin = {-64, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.UP_Seg uP_Seg3(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {-49, -43}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.PointMass LM3(m = 50, r_0(fixed = true, start = {0, -15.0, 0}), v_0(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {-49, -55}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    LMC_MBS.LO_Seg lO_Seg3(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {-49, -67}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(phi(fixed = true), useAxisFlange = true, w(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {-50, -26}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.DnB dnB4 annotation(
      Placement(visible = true, transformation(origin = {-26, 60}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Forces.WorldForce force3 annotation(
      Placement(visible = true, transformation(origin = {-4, 60}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.UP_Seg uP_Seg4(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {11, 73}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.PointMass LM4(m = 50, r_0(fixed = true, start = {0, -20.0, 0}), v_0(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {11, 61}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    LMC_MBS.LO_Seg lO_Seg4(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {11, 49}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4(phi(fixed = true), useAxisFlange = true, w(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {10, 90}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.DnB dnB5 annotation(
      Placement(visible = true, transformation(origin = {-26, 2}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Forces.WorldForce force4 annotation(
      Placement(visible = true, transformation(origin = {-4, 2}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.UP_Seg uP_Seg5(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {11, 15}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.PointMass LM5(m = 50, r_0(fixed = true, start = {0, -25.0, 0}), v_0(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {11, 3}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    LMC_MBS.LO_Seg lO_Seg5(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {11, -9}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5(phi(fixed = true), useAxisFlange = true, w(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {10, 32}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.DnB dnB6 annotation(
      Placement(visible = true, transformation(origin = {-26, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Forces.WorldForce force5 annotation(
      Placement(visible = true, transformation(origin = {-4, -56}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    LMC_MBS.UP_Seg uP_Seg6(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {11, -43}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.PointMass LM6(m = 50, r_0(fixed = true, start = {0, -30.0, 0}), v_0(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {11, -55}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    LMC_MBS.LO_Seg lO_Seg6(r = {0, -2.5, 0}) annotation(
      Placement(visible = true, transformation(origin = {11, -67}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6(phi(fixed = true), useAxisFlange = true, w(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {10, -26}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.PointMass Load(m = 100, r_0(fixed = true, start = {0, -32.5, 0}), v_0(fixed = true)) annotation(
      Placement(visible = true, transformation(origin = {11, -83}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Real Ef;
    Real X[7];
    Real Y[7];
  equation
    connect(LM3.frame_a, lO_Seg3.frame_a) annotation(
      Line(points = {{-48, -54}, {-48, -54}, {-48, -60}, {-48, -60}}, color = {95, 95, 95}));
    connect(lO_Seg6.frame_b, Load.frame_a) annotation(
      Line(points = {{12, -74}, {12, -74}, {12, -82}, {12, -82}}, color = {95, 95, 95}));
    connect(LM6.frame_a, lO_Seg6.frame_a) annotation(
      Line(points = {{12, -54}, {12, -54}, {12, -60}, {12, -60}}, color = {95, 95, 95}));
    connect(uP_Seg6.frame_b, LM6.frame_a) annotation(
      Line(points = {{12, -50}, {12, -50}, {12, -54}, {12, -54}}, color = {95, 95, 95}));
    connect(LM5.frame_a, lO_Seg5.frame_a) annotation(
      Line(points = {{12, 4}, {12, 4}, {12, -2}, {12, -2}}, color = {95, 95, 95}));
    connect(uP_Seg5.frame_b, LM5.frame_a) annotation(
      Line(points = {{12, 8}, {12, 8}, {12, 4}, {12, 4}}, color = {95, 95, 95}));
    connect(LM4.frame_a, lO_Seg4.frame_a) annotation(
      Line(points = {{12, 62}, {12, 62}, {12, 56}, {12, 56}}, color = {95, 95, 95}));
    connect(force5.frame_b, LM6.frame_a) annotation(
      Line(points = {{2, -56}, {12, -56}, {12, -54}, {12, -54}}, color = {95, 95, 95}));
    connect(dnB6.y1, force5.force) annotation(
      Line(points = {{-20, -56}, {-12, -56}, {-12, -56}, {-12, -56}}, color = {0, 0, 127}, thickness = 0.5));
    connect(lO_Seg6.y, dnB6.c2) annotation(
      Line(points = {{6, -72}, {-26, -72}, {-26, -62}, {-26, -62}}, color = {0, 0, 127}, thickness = 0.5));
    connect(uP_Seg6.y, dnB6.c1) annotation(
      Line(points = {{6, -38}, {-26, -38}, {-26, -50}, {-26, -50}}, color = {0, 0, 127}, thickness = 0.5));
    connect(revolute6.frame_b, uP_Seg6.frame_a) annotation(
      Line(points = {{16, -26}, {20, -26}, {20, -36}, {12, -36}, {12, -36}, {12, -36}}, color = {95, 95, 95}));
    connect(lO_Seg5.frame_b, revolute6.frame_a) annotation(
      Line(points = {{12, -16}, {0, -16}, {0, -26}, {4, -26}, {4, -26}, {4, -26}}, color = {95, 95, 95}));
    connect(force4.frame_b, LM5.frame_a) annotation(
      Line(points = {{2, 2}, {12, 2}, {12, 4}, {12, 4}}, color = {95, 95, 95}));
    connect(dnB5.y1, force4.force) annotation(
      Line(points = {{-20, 2}, {-12, 2}, {-12, 2}, {-12, 2}}, color = {0, 0, 127}, thickness = 0.5));
    connect(lO_Seg5.y, dnB5.c2) annotation(
      Line(points = {{6, -14}, {-26, -14}, {-26, -4}, {-26, -4}}, color = {0, 0, 127}, thickness = 0.5));
    connect(uP_Seg5.y, dnB5.c1) annotation(
      Line(points = {{6, 20}, {-26, 20}, {-26, 8}, {-26, 8}}, color = {0, 0, 127}, thickness = 0.5));
    connect(revolute5.frame_b, uP_Seg5.frame_a) annotation(
      Line(points = {{16, 32}, {20, 32}, {20, 22}, {10, 22}, {10, 22}, {12, 22}}, color = {95, 95, 95}));
    connect(lO_Seg4.frame_b, revolute5.frame_a) annotation(
      Line(points = {{12, 42}, {0, 42}, {0, 32}, {4, 32}, {4, 32}}, color = {95, 95, 95}));
    connect(lO_Seg4.y, dnB4.c2) annotation(
      Line(points = {{6, 44}, {-26, 44}, {-26, 54}, {-26, 54}}, color = {0, 0, 127}, thickness = 0.5));
    connect(uP_Seg4.frame_b, LM4.frame_a) annotation(
      Line(points = {{12, 66}, {12, 66}, {12, 62}, {12, 62}}, color = {95, 95, 95}));
    connect(force3.frame_b, LM4.frame_a) annotation(
      Line(points = {{2, 60}, {12, 60}, {12, 62}, {12, 62}}, color = {95, 95, 95}));
    connect(dnB4.y1, force3.force) annotation(
      Line(points = {{-20, 60}, {-12, 60}, {-12, 60}, {-12, 60}}, color = {0, 0, 127}, thickness = 0.5));
    connect(uP_Seg4.y, dnB4.c1) annotation(
      Line(points = {{6, 78}, {-26, 78}, {-26, 66}, {-26, 66}}, color = {0, 0, 127}, thickness = 0.5));
    connect(revolute4.frame_b, uP_Seg4.frame_a) annotation(
      Line(points = {{16, 90}, {20, 90}, {20, 80}, {12, 80}, {12, 80}, {12, 80}}, color = {95, 95, 95}));
    connect(lO_Seg3.frame_b, revolute4.frame_a) annotation(
      Line(points = {{-48, -74}, {-36, -74}, {-36, 90}, {4, 90}, {4, 90}}, color = {95, 95, 95}));
    connect(lO_Seg3.y, dnB3.c2) annotation(
      Line(points = {{-54, -72}, {-86, -72}, {-86, -62}, {-86, -62}}, color = {0, 0, 127}, thickness = 0.5));
    connect(uP_Seg3.frame_b, LM3.frame_a) annotation(
      Line(points = {{-48, -50}, {-48, -50}, {-48, -54}, {-48, -54}}, color = {95, 95, 95}));
    connect(force2.frame_b, LM3.frame_a) annotation(
      Line(points = {{-58, -56}, {-48, -56}, {-48, -54}, {-48, -54}}, color = {95, 95, 95}));
    connect(dnB3.y1, force2.force) annotation(
      Line(points = {{-80, -56}, {-72, -56}, {-72, -56}, {-72, -56}}, color = {0, 0, 127}, thickness = 0.5));
    connect(uP_Seg3.y, dnB3.c1) annotation(
      Line(points = {{-54, -38}, {-86, -38}, {-86, -50}, {-86, -50}, {-86, -50}}, color = {0, 0, 127}, thickness = 0.5));
    connect(revolute3.frame_b, uP_Seg3.frame_a) annotation(
      Line(points = {{-44, -26}, {-40, -26}, {-40, -36}, {-48, -36}, {-48, -36}, {-48, -36}}, color = {95, 95, 95}));
    connect(lO_Seg2.frame_b, revolute3.frame_a) annotation(
      Line(points = {{-48, -16}, {-60, -16}, {-60, -26}, {-56, -26}, {-56, -26}, {-56, -26}}, color = {95, 95, 95}));
    connect(lO_Seg2.y, dnB2.c2) annotation(
      Line(points = {{-54, -14}, {-86, -14}, {-86, -4}, {-86, -4}, {-86, -4}}, color = {0, 0, 127}, thickness = 0.5));
    connect(LM2.frame_a, lO_Seg2.frame_a) annotation(
      Line(points = {{-48, 4}, {-48, 4}, {-48, -2}, {-48, -2}}, color = {95, 95, 95}));
    connect(uP_Seg2.frame_b, LM2.frame_a) annotation(
      Line(points = {{-48, 8}, {-48, 8}, {-48, 4}, {-48, 4}}, color = {95, 95, 95}));
    connect(force1.frame_b, LM2.frame_a) annotation(
      Line(points = {{-58, 2}, {-48, 2}, {-48, 4}, {-48, 4}}, color = {95, 95, 95}));
    connect(dnB2.y1, force1.force) annotation(
      Line(points = {{-80, 2}, {-72, 2}, {-72, 2}, {-72, 2}}, color = {0, 0, 127}, thickness = 0.5));
    connect(uP_Seg2.y, dnB2.c1) annotation(
      Line(points = {{-54, 20}, {-86, 20}, {-86, 8}, {-86, 8}, {-86, 8}}, color = {0, 0, 127}, thickness = 0.5));
    connect(revolute2.frame_b, uP_Seg2.frame_a) annotation(
      Line(points = {{-44, 32}, {-40, 32}, {-40, 22}, {-50, 22}, {-50, 22}, {-48, 22}}, color = {95, 95, 95}));
    connect(lO_Seg1.frame_b, revolute2.frame_a) annotation(
      Line(points = {{-48, 42}, {-60, 42}, {-60, 32}, {-56, 32}, {-56, 32}, {-56, 32}}, color = {95, 95, 95}));
    connect(LM1.frame_a, lO_Seg1.frame_a) annotation(
      Line(points = {{-48, 62}, {-48, 62}, {-48, 56}, {-48, 56}}, color = {95, 95, 95}));
    connect(uP_Seg1.frame_b, LM1.frame_a) annotation(
      Line(points = {{-48, 66}, {-48, 66}, {-48, 62}, {-48, 62}}, color = {95, 95, 95}));
    connect(force.frame_b, LM1.frame_a) annotation(
      Line(points = {{-58, 60}, {-50, 60}, {-50, 62}, {-48, 62}}, color = {95, 95, 95}));
    connect(dnB1.y1, force.force) annotation(
      Line(points = {{-80, 60}, {-72, 60}, {-72, 60}, {-72, 60}}, color = {0, 0, 127}, thickness = 0.5));
    connect(lO_Seg1.y, dnB1.c2) annotation(
      Line(points = {{-54, 44}, {-86, 44}, {-86, 54}, {-86, 54}}, color = {0, 0, 127}, thickness = 0.5));
    connect(uP_Seg1.y, dnB1.c1) annotation(
      Line(points = {{-54, 78}, {-86, 78}, {-86, 66}, {-86, 66}}, color = {0, 0, 127}, thickness = 0.5));
    connect(revolute1.frame_b, uP_Seg1.frame_a) annotation(
      Line(points = {{-44, 90}, {-40, 90}, {-40, 80}, {-48, 80}, {-48, 80}}, color = {95, 95, 95}));
    connect(fixed1.frame_b, revolute1.frame_a) annotation(
      Line(points = {{-64, 90}, {-56, 90}, {-56, 90}, {-56, 90}}, color = {95, 95, 95}));
    Ef = sqrt(fixed1.frame_b.f[1] ^ 2 + fixed1.frame_b.f[2] ^ 2);
    X[1]=uP_Seg1.frame_a.r_0[1];
    X[2]=uP_Seg2.frame_a.r_0[1];
    X[3]=uP_Seg3.frame_a.r_0[1];
    X[4]=uP_Seg4.frame_a.r_0[1];
    X[5]=uP_Seg5.frame_a.r_0[1];
    X[6]=uP_Seg6.frame_a.r_0[1];
    X[7]=lO_Seg6.frame_b.r_0[1];
    Y[1]=uP_Seg1.frame_a.r_0[2];
    Y[2]=uP_Seg2.frame_a.r_0[2];
    Y[3]=uP_Seg3.frame_a.r_0[2];
    Y[4]=uP_Seg4.frame_a.r_0[2];
    Y[5]=uP_Seg5.frame_a.r_0[2];
    Y[6]=uP_Seg6.frame_a.r_0[2];
    Y[7]=lO_Seg6.frame_b.r_0[2];
    annotation(
      experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.1));
  end SubseaLoad;


  annotation(
    Diagram(coordinateSystem(initialScale = 0.05)),
    uses(Modelica(version = "3.2.2")));
end LMC_MBS;
