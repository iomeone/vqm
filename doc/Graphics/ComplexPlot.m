(* :Title:   ComplexPlot *)

(* :Name:    Graphics`ComplexPlot` *)

(* :Author:  Bernd Thaller,
             Institute of Mathematics,
             University of Graz,
             A-8010 Graz
             bernd.thaller@kfunigraz.ac.at
*)

(* :Date:    Dec 27, 1999 *)

(* :Summary:
    Generate two-dimensional ColorDensityGraphics
    and three-dimensional SurfaceGraphics from
    complex-valued functions with a one-to-one
    color mapping for complex numbers.
*)

(* :Package Version:        1.8 *)

(* :Mathematica Version:    3.0 and 4.01 *)

(* :Copyright: Copyright 1998, 1999 Bernd Thaller *)

(* :Keywords:
    DensityPlot, Complex Function, Wavefunction
*)

(* :Limitations:
    The Background option does not work properly.
    Compile errors with Mathematica Versions
    prior to 3.0.
    'ComplexContourPlot' is just a shorthand for combining
    a ComplexDensityPlot with a ContourPlot of Abs[f].
    Only 'ListComplexSurfacePlot' is available in
    this version (functions have to be converted to
    a table of values first, see the examples),
    ListComplexSurfacePlot has no Mesh-option.
*)

(* :Acknowledgement:
    Thanks to Manfred Liebmann for some
    performance-improving suggestions.
*)

(*-----------------------------------*)
BeginPackage[
    "Graphics`ComplexPlot`",    (* package Context *)
    "Graphics`FilledPlot`",     (* draw axes in front *)
    "Utilities`FilterOptions`",
    "Global`"]                  (* prevent shadowing *)
(*-----------------------------------*)

Off[General::spell1,General::spell];

Clear[$ComplexToColorMap];

$MaxAbsValue=1;

Unprotect[ComplexPlot3D,ListComplexPlot3D,
ListComplexSurfacePlot,ComplexDensityPlot,
ListComplexDensityPlot,ComplexContourPlot,
ListComplexContourPlot,ColorArrayPlot,ColorDensityGraphics,
ComplexToColor,ComplexToColorMap,ValueRange,LightnessRange,
SphereRadius,Highlights,ValueChecking,ScaledValues]

ComplexPlot::usage =
"This package provides commands for visualizing complex-valued
functions by generating two-dimensional ColorDensityGraphics,
ContourGraphics and three dimensional SurfaceGraphics
of complex-valued functions with a color code for complex
numbers.";

ComplexPlot3D::usage =
"ComplexPlot3D[f[x,y],{x,xmin,xmax},{y,ymin,ymax},opts]
generates a surface plot of a complex-valued function f of
two variables. The height of the surface is given by the
absolute value, the color is determined by the complex value
of f according to the option ComplexToColorMap.";

ListComplexPlot3D::usage =
"ListComplexPlot3D[array,opts] generates a SurfaceGraphics
from a two-dimensional array of complex numbers. The height
of the surface is given by the absolute value, the color is
determined by the complex value of f according to the option
ComplexToColorMap.";

ListComplexSurfacePlot::usage =
"Similar to ListComplexPlot3D, but with a 'real surface look'.
The option Highlights->On (default) lets the surface appear
glossy.";

ComplexDensityPlot::usage =
"ComplexDensityPlot[f[x,y],{x,xmin,xmax},{y,ymin,ymax},opts]
generates a ColorDensityGraphics of a complex-valued function f
of two variables x and y. It is similar to DensityPlot. The
complex value f[x,y] is mapped one-to-one to a color. The color
map is given by the option ComplexToColorMap. The default
$ComplexToColorMap determines the hue of the color from Arg[z],
and the lightness from Abs[z].";

ListComplexDensityPlot::usage =
"ListComplexDensityPlot[array,opts] gives a ColorDensityGraphics
of a two-dimensional array of complex numbers. It is similar to
ListDensityPlot. Each complex number is mapped one-to-one to a
color. The color map is determinded by the option
ComplexToColorMap. The default $ComplexToColorMap determines
the hue of the color from Arg[z], and the lightness from
Abs[z].";

ComplexContourPlot::usage =
"ComplexContourPlot[f[x,y],{x,xmin,xmax},{y,ymin,ymax},opts]
visualizes a complex-valued function f of two variables
x and y. ComplexContourPlot combines a ColorDensityGraphics
with a ContourGraphics of the absolute value of f.";

ListComplexContourPlot::usage =
"ListComplexContourPlot[array,opts] generates a
ColorDensityGraphics of a two-dimensional array of complex
numbers and combines it with a ContourGraphics of Abs[array].";

ColorArrayPlot::usage =
"ColorArrayPlot[list, opts] makes a 2D raster plot with colors
given by list. Here list is a two-dimensional array of color
directives.";

ComplexToColor::usage =
"ComplexToColor[z,opts] associates a color to a complex number z.
The color map is given by the option ComplexToColorMap.";

ComplexToColorMap::usage =
"ComplexToColorMap is an option for ComplexToColor,
ComplexDensityPlot, and ComplexPlot3D. The default setting is
ComplexToColorMap->$ComplexToColorMap, which determines
the color representing a complex number via a stereographic
projection from the complex plane onto the surface of the
color manifold in the HLS system";

ValueRange::usage =
"ValueRange is an option for ComplexToColor, ComplexDensityPlot,
and ComplexPlot3D which is used to determine parameters for the
default color map $ComplexToColorMap. ValueRange->{rmin,rmax},
where rmin < rmax are positive real numbers, causes every z with
Abs[z]<rmin (resp. Abs[z]>rmax) to be colored with minimal (resp.
maximal) lightness. If rmax < rmin, then Abs[z] < rmax has maximal
lightness, Abs[z] > rmin has minimal lightness. The default
Automatic corresponds to ValueRange->{0,Infinity}. The setting
for minimal/maximal lightness is determined by the
LightnessRange option";

LightnessRange::usage =
"LightnessRange is an option for ComplexToColor, ComplexDensityPlot,
and ComplexPlot3D. LightnessRange->{lmin,lmax}, where lmin and lmax
are real numbers in the interval [0,1], sets the minimal lightness
to lmin and the maximal lightness to lmax. Default is
LightnessRange->{0,1}. The setting for this option is used to
determine parameters for the default color map $ComplexToColorMap.";

SphereRadius::usage =
"SphereRadius is an option for ComplexToColor, ComplexDensityPlot,
and ComplexPlot3D. SphereRadius->R sets the radius of the sphere
used by $ComplexToColorMap to R. $ComplexToColorMap uses a
stereographic projection onto the surface of the color manifold
to determine the color of a complex number. Setting the radius
to R causes complex numbers with Abs[z]=R to be drawn at
lightness 1/2 (which corresponds to maximal brightness and
saturation in the HSB system).";

ValueChecking::usage =
"ValueChecking is an option for ComplexDensityPlot and
ComplexContourPlot. It has no effect on the Plot3D commands.
The complex values 0 and ComplexInfinity have no well defined
color, in particular when LightnessRange is {lmin,lmax} instead
of {0,1}. The default ValueChecking->On plots 0 as
GrayLevel[lmin], ComplexInfinity as GrayLevel[lmax] and
Indeterminate values in an intermediate gray level.
For slightly better performance in the case that there are
no exceptional points, use the setting ValueChecking->Off.";

ScaledValues::usage = "Option for ComplexDensityPlot etc.
If ScaledValues is set to True, then all values are scaled
so that the maximal absolute value is $MaxAbsValue (=1 by default).";

Highlights::usage =
"Highlights is an option for ListComplexSurfacePlot.";

ColorDensityGraphics::usage =
"ColorDensityGraphics[absarray,colorarray,{opts}] is a
representation of a two-dimensional plot of an array of complex
numbers. It can be converted to SurfaceGraphics, ContourGraphics,
DensityGraphics and Graphics objects.";

$ComplexToColorMap::usage =
"$ComplexToColorMap[r,phi,{parameters}],
associates a color to a complex number. The complex number is
given in polar form, r=Abs[z], phi=Arg[z]. The Hue of the color
is given by phi/2/Pi, the lightness is determined from r. The
color map can be described as a stereographic projection from
the complex plane onto the surface of the color manifold in the
Hue-Lightness-Saturation system. The optional parameters are
{R,bmin,bmax,lmin,lmax} specifying the radius R of the sphere, the
values bmin, bmax depend on the values of r which are to be drawn
with minimal resp. maximal lightness lmin resp. lmax.";

(*-----------------------------------*)
Begin["`Private`"]
(*-----------------------------------*)

myoptions =
    {ComplexToColorMap->Automatic,
        LightnessRange->Automatic,
        SphereRadius->1., ValueRange->Automatic, ScaledValues->False,
        ValueChecking->On, Highlights->On};


(* --------- ComplexPlot3D --------- *)

Options[ComplexPlot3D] =
    Join[myoptions,Options[Plot3D]];

SetAttributes[ComplexPlot3D,HoldAll]

ComplexPlot3D[func_,{x_Symbol,xmin_,xmax_},
        {y_Symbol,ymin_,ymax_}, opts:((_Rule | _RuleDelayed)...)] := 
Module[{fli,cls1,cls,
        colormap=ComplexToColorMap/.{opts}/.myoptions},
       If[colormap === Automatic || colormap === $ComplexToColorMap,
          fli = clight[processOptions[opts]];
          cls1 = Hue[mod1[#2], sat[#1], bri[#1]]&;
          cls = cls1[fli[#1],#2]&,
       (*else*) cls = colormap];
       Plot3D[{Abs[func], cls[Abs[func],Arg[func]]},
              {x,xmin,xmax}, {y,ymin,ymax},
              Evaluate[FilterOptions[Plot3D,opts]]]
]/;And @@ NumberQ /@ N[{xmin,xmax,ymin,ymax}];


(* ------- ListComplexPlot3D ------- *)

Options[ListComplexPlot3D] =
    Join[myoptions, Options[ListPlot3D]];
    
ListComplexPlot3D[array_, opts:((_Rule | _RuleDelayed)...)] := 
Module[{ab1=Abs[array],absarr,
        argarr=N[Arg[Drop[#,-1]& /@ Drop[array,-1]]],args,
        abss,ligarr,hues,stns,brts,colors,
        colormap=ComplexToColorMap/.{opts}/.myoptions,
        scaledvl=ScaledValues/.{opts}/.myoptions},
       If[scaledvl === True,absarr = $MaxAbsValue ab1/Max[ab1],absarr = ab1];
       Remove[ab1];
       abss=Drop[#,-1]& /@ Drop[absarr,-1];
       If[colormap===Automatic || colormap === $ComplexToColorMap,
             args = Map[Replace[#,a_/;!NumberQ[a] ->0]&,argarr,{2}];
             Remove[argarr];
             ligarr = Map[clight[processOptions[opts]], abss, {2}];
             Remove[abss];
             hues = Map[mod1, args, {2}];
	         Remove[args];
	         stns = Map[sat, ligarr, {2}];
	         brts = Map[bri, ligarr, {2}];
	         Remove[ligarr];
	         colors = MapThread[Hue, {hues,stns,brts}, 2],
	         Remove[hues,stns,brts],
       (*else*) colors = MapThread[colormap,{abss,args},2];
             Remove[abss,args]];
       ListPlot3D[ absarr, colors,
            Evaluate[FilterOptions[ListPlot3D,opts]]
       ]
]/; MatrixQ[array]

(* ------- ListComplexSurfacePlot ------- *)

Options[ListComplexSurfacePlot] =
    Join[myoptions, Options[Graphics3D]];

ListComplexSurfacePlot[array_, opts:((_Rule | _RuleDelayed)...)] := 
Module[{abss = Abs[array], polys,
        args = Map[ Replace[#,a_/;!NumberQ[a] ->0]&,
                   N[Arg[Drop[#,-1]& /@ Drop[array,-1]]], {2}],
        clgs,hues,stns,brts,colors,
        hlite = Highlights /.{opts}/.Options[ListComplexSurfacePlot]},
       polys =
         Graphics3D[
           ListPlot3D[abss,DisplayFunction->Identity]
         ][[1]];
       clgs = Map[ clight[processOptions[opts]],
                   Drop[#,-1]& /@ Drop[abss,-1], {2}];
       Remove[abss];
       hues = Map[mod1, args, {2}];
	   Remove[args];
	   stns = Map[sat, clgs, {2}];
	   brts = Map[bri, clgs, {2}];
	   Remove[clgs];
	   If[hlite === On,
	   colors =
	          Flatten[
	            MapThread[
	             SurfaceColor[Hue[#1,#2,#3],GrayLevel[1.],15]& ,
	            {hues,stns,brts}, 2]],
	   colors =
	          Flatten[
	            MapThread[
	             SurfaceColor[Hue[#1,#2,#3]]& ,
	            {hues,stns,brts}, 2]]
	   ];
	   Remove[hues,stns,brts];
       Show[
         Graphics3D[{EdgeForm[],{colors,polys}//Transpose}],
            Evaluate[FilterOptions[Graphics3D,opts]],
            LightSources->{{{1,1,1},GrayLevel[1]},
                           {{1,0,1},GrayLevel[.4]},
                           {{1,-1,1},GrayLevel[.3]}},
            BoxRatios->{1,1,.4},
            Axes->True
       ]
]/; MatrixQ[array]


(* ------- ComplexDensityPlot ------- *)

Options[ComplexDensityPlot] =
    Join[myoptions, Options[DensityPlot],
    {AxesFront->True}];

SetAttributes[ComplexDensityPlot,HoldAll]

ComplexDensityPlot[func_,
        {x_Symbol,xmin_,xmax_}, {y_Symbol,ymin_,ymax_},
        opts:((_Rule | _RuleDelayed)...)] :=
    Module[{pl, comp, ops, dx, f=Function[{x,y},func],
            bmax=N[{xmax,ymax}],bmin=N[{xmin,ymin}], gc, array},
        {pl,comp} = {PlotPoints,Compiled}/.{opts}/.Options[DensityPlot];
        pl = testplotpoints[pl];
        ops = JoinOptions[MeshRange->{{xmin,xmax},
                {ymin,ymax}},opts];
        dx  = (bmax-bmin)/(pl-1.);
        If[comp == True,
           gc = fcomp[f,bmin,dx];
           array =
              Table[gc[j,i],{j,0,pl[[2]]-1},{i,0,pl[[1]]-1}],
           array = Table[f[x,y],
                {y, N[ymin], ymax, dx[[2]]},
                {x, N[xmin], xmax, dx[[1]]}] ];
        ListComplexDensityPlot[array,ops]
    ]/; And @@ NumberQ /@ N[{xmin,xmax,ymin,ymax}]


fcomp[f_, {xmin_, ymin_}, {dx_, dy_}] := 
	Compile @@ {{{j,_Integer},{i,_Integer}},
	Simplify[f[(xmin+i*dx),(ymin+j*dy)]]}


(* ---- ListComplexDensityPlot ---- *)

Options[ListComplexDensityPlot] =
    Join[myoptions, Options[ListDensityPlot],
    {AxesFront->True}];
    
ListComplexDensityPlot[array_, opts:((_Rule | _RuleDelayed)...)] :=
    Module[{arr=N[array],ab1,absarr,argarr,ligarr,hues,stns,brts,
            colors,params=processOptions[opts],
            checkopts=
              {ComplexToColorMap,ValueChecking}/.{opts}/.myoptions,
            scaledvl=ScaledValues/.{opts}/.myoptions},
       ab1 = Abs[arr];
       If[scaledvl === True, absarr = $MaxAbsValue ab1/Max[ab1], absarr = ab1];
       Remove[ab1];
       argarr = Map[Replace[#,a_/;!NumberQ[a] ->0]&,N[Arg[arr]],{2}];
       If[checkopts[[1]] === Automatic ||
          checkopts[[1]] === $ComplexToColorMap,
             Off[CompiledFunction::"cfr"];
             ligarr = Map[clight[params], absarr, {2}];
             hues   = Map[mod1, argarr, {2}];
	         Remove[argarr];
	         stns   = Map[sat, ligarr, {2}];
	         brts   = Map[bri, ligarr, {2}];
	         Remove[ligarr];
             On[CompiledFunction::"cfr"];
	         colors = MapThread[Hue, {hues,stns,brts}, 2];
	         Remove[hues,stns,brts];
	         If[checkopts[[2]]===On || checkopts[[2]] === True,
	         colors = checkvalues[arr,colors,params]; Remove[arr],
	         Remove[arr]],
          (*else*) colors =
            MapThread[checkopts[[1]][#1,#2]&,{absarr,argarr},2];
            Remove[argarr]];
       Show[ColorDensityGraphics[ absarr, colors,
        {JoinOptions[{opts},Options[ColorDensityGraphics]]}]]
    ]/; MatrixQ[N[array]]


checkvalues[arr_,colors_,params_]:=
   Module[{colors1=colors,bmax=params[[3]],
           lmin=params[[4]],lmax=params[[5]]},
      {li,la}=If[bmax<=0,{lmax,lmin},{lmin,lmax}];
      colors1=ReplacePart[colors1, Hue[0,0,li],
                      Position[arr,z_/;z==0]];
      colors1=ReplacePart[colors1, Hue[0,0,.5 (lmin+lmax)],
                      Position[arr,z_/;z==Indeterminate]];
      ReplacePart[colors1,Hue[0,0,la],
                      Position[arr,z_/;z==ComplexInfinity]]                   

   ];


(* ------ ColorDensityGraphics ------ *)

Options[ColorDensityGraphics] =
    Join[Options[ListDensityPlot],{AxesFront->True}];

Format[_ColorDensityGraphics] = "-ColorDensityGraphics-";

ColorDensityGraphics/:
Show[ColorDensityGraphics[abs_,colors_,opts___],opts2___] :=
    Module[{ops=FilterOptions[
            ColorArrayPlot,JoinOptions[{opts2},opts]]},
        ColorArrayPlot[colors,ops];
        ColorDensityGraphics[abs,colors,{ops}]
    ]

ColorDensityGraphics[SurfaceGraphics[abs_,colors_,opts___],opts2___] :=
    ColorDensityGraphics[abs,colors,
    FilterOptions[ColorDensityGraphics,JoinOptions[{opts2},opts]]]

ColorDensityGraphics/:
SurfaceGraphics[ColorDensityGraphics[abs_,colors_,opts___],opts2___] :=
    Module[{colorarray},
        If[Dimensions[abs]==Dimensions[colors],
            colorarray = Drop[#,-1]& /@ Drop[colors,-1],
            colorarray = colors];
        SurfaceGraphics[abs, colorarray,
            FilterOptions[SurfaceGraphics,JoinOptions[{opts2},opts]]]
    ]

ColorDensityGraphics/:
DensityGraphics[ColorDensityGraphics[abs_,colors_,opts___],opts2___] :=
    DensityGraphics[abs,
        FilterOptions[DensityGraphics,JoinOptions[{opts2},opts]]]

ColorDensityGraphics/:
ContourGraphics[ColorDensityGraphics[abs_,colors_,opts___],opts2___] :=
    ContourGraphics[abs,
        FilterOptions[ContourGraphics,JoinOptions[{opts2},opts]]]

ColorDensityGraphics/:
Graphics[ColorDensityGraphics[abs_,colors_,opts___]] :=
    ColorArrayGraphics[colors,opts]


(* -------- ColorArrayPlot --------- *)

Options[ColorArrayPlot] = Options[ColorDensityGraphics];

ColorArrayPlot[colors_,opts:((_Rule | _RuleDelayed)...)] :=
    Show[ColorArrayGraphics[colors,opts]]

ColorArrayGraphics[colors_,opts___] :=
    Module[{mshran,msh,mshsty,dims=Dimensions[colors],
           ops=JoinOptions[Flatten[{opts}],Options[ColorArrayPlot]]},
        {mshran,msh,mshsty}={MeshRange,Mesh,MeshStyle}/.{ops};
        If[ mshran === Automatic,
            mshran = { {0,dims[[2]]}, {0,dims[[1]]} } ];
        ops=FilterOptions[Graphics,ops];
        If[ mshsty === Automatic, mshsty = GrayLevel[0] ];
        If[msh === True || msh === Automatic,
            Graphics[
                {RasterArray[colors,mshran//Transpose],
                Flatten[ List @@
                meshlines[mshran,{dims[[2]],dims[[1]]},mshsty]]},
                {ops}],
            Graphics[RasterArray[colors,mshran//Transpose],{ops}]
        ]
   ]


(* ------- ComplexContourPlot ------- *)

Options[ComplexContourPlot]=
    Join[myoptions,Options[ContourPlot]];

SetAttributes[ComplexContourPlot,HoldAll];

ComplexContourPlot[func_,
        {x_Symbol,xmin_,xmax_},{y_Symbol,ymin_,ymax_},
        opts:((_Rule | _RuleDelayed)...)]:=
Module[{gr1,gr2},
    gr1=ComplexDensityPlot[func,{x,xmin,xmax},{y,ymin,ymax},
        Mesh->False,PlotRange->All,DisplayFunction->Identity,opts];
    gr2=ContourGraphics[gr1,ContourShading->False,
        FilterOptions[ContourGraphics,opts]];
    Show[Graphics[gr1],gr2,PlotRange->All,FilterOptions[Graphics,opts],
        DisplayFunction->$DisplayFunction]
]/; And @@ NumberQ /@ N[{xmin,xmax,ymin,ymax}]


(* ----- ListComplexContourPlot ----- *)

Options[ListComplexContourPlot]=
    Join[myoptions,Options[ListContourPlot]];

ListComplexContourPlot[array_,opts:((_Rule | _RuleDelayed)...)]:=
Module[{gr1,gr2,dims=Dimensions[array]},
    gr1=ListComplexDensityPlot[array,
        Mesh->False,PlotRange->All,DisplayFunction->Identity,opts];
    gr2=ContourGraphics[gr1,ContourShading->False,
        MeshRange->{ {0,dims[[2]]}, {0,dims[[1]]} },
        FilterOptions[ContourGraphics,opts]];
    Show[Graphics[gr1],gr2,PlotRange->All,FilterOptions[Graphics,opts],
        DisplayFunction->$DisplayFunction]
]/; MatrixQ[N[array]]


(* --------- ComplexToColor -------- *)

Options[ComplexToColor] = myoptions;

SetAttributes[ComplexToColor,Listable]

ComplexToColor[z_,opts:((_Rule | _RuleDelayed)...)] :=
    Module[{colormap=ComplexToColorMap/.{opts}/.myoptions},
        If[colormap===Automatic ||
        colormap === $ComplexToColorMap,
        $ComplexToColorMap[Abs[z],Arg[z],processOptions[opts]],
        colormap[Abs[z],Arg[z]]]
    ]


(* ------ $ComplexToColorMap ------- *)

$ComplexToColorMap[r_,arg_,
        {R_:1.,bmin_:0.,bmax_:1.,lmin_:0.,lmax_:1.}] :=
    Hue[0, 0, If[bmax<0,lmax,lmin]]/;r==0

$ComplexToColorMap[Infinity,arg_,
        {R_:1.,bmin_:0.,bmax_:1.,lmin_:0.,lmax_:1.}] :=
    Hue[mod1[arg], sat[#], bri[#]]&[If[bmax<0,lmin,lmax]]/;NumberQ[N[arg]]

$ComplexToColorMap[Infinity,arg_,
        {R_:1.,bmin_:0.,bmax_:1.,lmin_:0.,lmax_:1.}] :=
    Hue[0, 0,If[bmax<0,lmin,lmax]]/;!NumberQ[N[arg]]

$ComplexToColorMap[r_,arg_,
        {R_:1.,bmin_:0.,bmax_:1.,lmin_:0.,lmax_:1.}] :=
    Hue[mod1[arg], sat[#], bri[#]
    ]&[light[r,R,bmin,bmax,lmin,lmax]]/;NumberQ[N[r]]/;NumberQ[N[r]]

$ComplexToColorMap[r_,arg_,
        {R_:1.,bmin_:0.,bmax_:1.,lmin_:0.,lmax_:1.}] :=
    Hue[0, 0,(lmin+lmax)/2]/;!NumberQ[N[r]]

light = Compile[{r,R,bmin,bmax,lmin,lmax},
           Min[{1,Max[{0.,2/Pi ArcTan[r/R] (bmax-bmin) + bmin}]}](
            lmax-lmin)+lmin ];

clight[{R_:1.,min_:0.,max_:1.,lmin_:0.,lmax_:1.}]:=
    Compile @@
     {{r},Simplify[Min[{1,Max[{0.,2/Pi ArcTan[r/R] (max-min) + min}]}](
            lmax-lmin)+lmin]};

mod1 = Compile[{arg},
	If[arg < 0.0, 1.0 + arg/(2 Pi), arg/(2 Pi)]];

sat   = Compile[{li}, If[li<=.5, 1, 2-2 li]];

bri   = Compile[{li}, If[li<=.5, 2 li, 1]];


(* determine params for $ComplexToColorMap from options *)

processOptions[opts___] :=
    Module[ {rad,vran,lran,bds},
        {rad,vran,lran} =
            {SphereRadius,ValueRange,LightnessRange}
                /.{opts}/.myoptions;
        rad = testrad[rad];
        bds = borders[ testvran[vran]/rad];
        Flatten[{rad,bds,testlran[lran]}]
    ]

borders[{0.,Infinity}] = {0.,1.};
borders[{Infinity,0.}] = {1.,0.};
borders[vran_] :=
    ({#[[1]],#[[1]]+N[Pi]})/(#[[1]]-#[[2]])&[-2 ArcTan[vran]]


(* -------- testing options --------- *)

testplotpoints[pl_] :=
    Module[{},
    	If[ IntegerQ[pl] && Positive[pl], {pl,pl}, pl]
    	] 

testrad[rad_]:=
    Module[{},
        If[ !NumberQ[N[rad]] || !Positive[N[rad]],
            Message[ComplexPlot::badradius]; Return[1.],
            Return[N[rad]]
        ]
    ]

testvran[vran_]:=
    Module[{},
        If[vran===Automatic || vran=={0,Infinity},
            Return[{0.,Infinity}]];
        If[vran=={Infinity,0},Return[{Infinity,0.}]];
        If[ !VectorQ[vran] || Length[vran] != 2 || 
            Negative[N[vran[[1]]]] ||
            Negative[N[vran[[2]]]] ||
            N[vran[[1]]]==N[vran[[2]]],
            Message[ComplexPlot::badvrange]; {0.,Infinity},
            N[vran] ]
    ]

testlran[lran_]:=
    Module[{},
        If[lran===Automatic,Return[{0.,1.}]];
        If[ !VectorQ[lran]  || Length[lran] != 2 || 
            N[lran[[1]]]<0. || N[lran[[2]]]>1.   ||
            N[lran[[1]]]>1. || N[lran[[2]]]<0. ,
            Message[ComplexPlot::badlrange];{0.,1.},
            N[lran]]
    ]


(* ------ dealing with options ------ *)

(* joining two lists of rules with list1 having precedence
over list2: *)

JoinOptions[list1_,list2___]:=
    Module[{namelist=First /@ Flatten[{list1}]},
        Sequence @@
        Sort[Join[Flatten[{list1}],
            Select[Flatten[{list2}],
                !MemberQ[namelist,First[#]]&]]
        ]
    ]


(* ------ generate mesh lines ------ *)

meshlines[m_,dms_,style_] :=
    Graphics[
        Flatten[{{style},
        Table[Line[{{i,m[[2,1]]},{i,m[[2,2]]}}],
            {i,m[[1,1]],m[[1,2]],(m[[1,2]]-m[[1,1]])/dms[[1]]}],
        Table[Line[{{m[[1,1]],i},{m[[1,2]],i}}],
            {i,m[[2,1]],m[[2,2]],(m[[2,2]]-m[[2,1]])/dms[[2]]}]
        }]
    ]


(* ----------- messages ------------ *)

ComplexPlot::badvrange =
"ValueRange should be Automatic or a list of two different
positive numbers (including 0 and Infinity).
Using Automatic instead.";

ComplexPlot::badlrange =
"LightnessRange should be Automatic or a list of two numbers
between 0 and 1. Using Automatic instead.";

ComplexPlot::badradius =
"SphereRadius must be a positive real number.
Using SphereRadius->1. instead.";

(*-----------------------------------*)
End[]       (* end `Private` Context *)
(*-----------------------------------*)

Protect[ComplexPlot3D,ListComplexPlot3D,
ListComplexSurfacePlot,ComplexDensityPlot,
ListComplexDensityPlot,ComplexContourPlot,
ListComplexContourPlot,ColorArrayPlot,ColorDensityGraphics,
ComplexToColor,ComplexToColorMap,ValueRange,LightnessRange,
SphereRadius,Highlights,ValueChecking,ScaledValues]

On[General::spell1,General::spell];

(*-----------------------------------*)
EndPackage[]  (* end package Context *)
(*-----------------------------------*)

(* :Examples:

ComplexDensityPlot[x+I y, {x,-4,4}, {y,-4,4}]

gr1=ComplexDensityPlot[Sin[x + I y],
        {x,-Pi,Pi}, {y,-1,2.5},
        Mesh->False, PlotPoints->30];

Show[SurfaceGraphics[gr1],
    AspectRatio->Automatic,Axes->True,Mesh->True]

ComplexContourPlot[x + I y, {x,-4,4}, {y,-4,4},
    ValueRange->{2,5},
    PlotRange->{0,3},Contours->5,
    ContourStyle->GrayLevel[0.5]]

ComplexContourPlot[Zeta[x + I y],
    {x,-.7,2.5}, {y,-2,42},
    PlotPoints->{10,50}, PlotRange->{-1.5,8},
    Contours->8, AspectRatio->5]

myfunc[{y_,m_,c_}] := Hue[-2 y,1-m,c];
collist=Table[myfunc[Mod[{x,x+y,x-y},1]],
        {y,0,1,1/20.}, {x,0,1,1/20.}];
ColorArrayPlot[collist,Mesh->False]

lis = Table[N[Tan[x + I y]],
        {y,-1.,1,.1},{x,-N[Pi],Pi,N[Pi]/10}];
ListComplexDensityPlot[lis,
    MeshRange->{{-Pi,Pi},{-1,1}}]

ComplexPlot3D[2 Exp[-x^2 - y^2 - 3 I x] Sin[x + I y],
    {x,-1,1},{y,-1,1}]

f[x_,y_]:=
   Which[Abs[x+I y] > 1.1, Indeterminate,
		Abs[x+I y] > 0.8, 0,
		Abs[x+I y] > 0.5, ComplexInfinity,
		True, DirectedInfinity[x+I y]];
Off[CompiledFunction::"cfr",Graphics::realu,Arg::indet ];
ComplexDensityPlot[f[x,y],{x,-1,1},{y,-1,1},
    Compiled->False, LightnessRange->{.2,.8}];
ComplexDensityPlot[f[x,y],{x,-1,1},{y,-1,1},
    ValueChecking->Off,
    Compiled->False, LightnessRange->{.2,.8}];
On[CompiledFunction::cfr,Graphics::realu,Arg::indet ];

arr = Table[Exp[I 3 (x+y) - x^2-y^2] +
            Exp[-I 2 (x+y) - (x-.5)^2-(y-.5)^2],
               {y,-3,3,.1},{x,-3,3,.1}];
ListComplexSurfacePlot[arr,
    PlotRange->All, SphereRadius->0.5,
    LightnessRange->{0.2,1.}];               

*)
