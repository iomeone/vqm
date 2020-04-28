(* :Title:	Arg Color Plot *)

(* :Name:	Graphics`ArgColorPlot` *)

(* :Author:	Bernd Thaller,
			Institute of Mathematics,
			University of Graz,
			A-8010 Graz
			bernd.thaller@kfunigraz.ac.at
*)

(* :Summary: 
	Plot the absolute value Abs[f[x]] of a complex-valued
	function f depending on a real variable x and
	fill the space between the plotted function
	and the x-axis with a color corresponding to the argument
	Arg[f[x]]. The saturation and brightness
	of the colors can be set using the options Saturation
	and Brightness.
*)

(* :Package Version:		1.3 *)

(* :Mathematica Version:	2.2.2 and 3.0 *)

(* :Copyright: Copyright 1997, Bernd Thaller *)

(* :Keywords:
	graphics, Plot, area under the curve, ArgColors
*)

(* :Limitations:
    Only one complex function can be plotted at a time, but
    a plot of a complex function can be combined with a plot
    of a real-valued function.
    For speed reasons, the complex function is compiled. The
    compiler may sometimes display error messages, but plotting
    still works.
*)
    
(* :Examples:
    Try the following:

    ArgColorPlot[Exp[-I 6 x - x^2/2], {x,-4,4}];

    ArgColorPlot[Exp[-I 6 x - x^2/2], {x,-4,4},
           PlotStyle -> {Thickness[0.025],GrayLevel[0.5]},
           Saturation -> .5, Brightness->.7,
           PlotRange -> {-.2,1.2}, Frame -> True,
           Axes->{True,False}]
    
    mytab = Table[Cos[x] + I Sin[x], {x,-4,4,.2}];
    ListArgColorPlot[mytab];
    
    CombinedPlot[{Gamma[I+x],Sin[Pi x]}, {x,-4,3}];

    mytab = Table[Sin[Pi x], {x,-4,3,.1}];
    ListCombinedPlot[{mytab,Sin[Pi x]}, {x,-4,3}];
*)
    
BeginPackage[
    "Graphics`ArgColorPlot`",
    "Graphics`FilledPlot`",     (* draw axes in front *)
    "Utilities`FilterOptions`",
    "Graphics`Graphics`",
    "Global`"]                  (* prevent shadowing *)

Unprotect[ArgColorPlot,
    ListArgColorPlot, CombinedPlot,
    ListCombinedPlot, NiceTicks];

Clear[ArgColorPlot,
    ListArgColorPlot, CombinedPlot,
    ListCombinedPlot, NiceTicks];
    
ArgColorPlot::usage = "ArgColorPlot[f[x],{x,x0,x1},opts] is used
like the usual Plot command. It gives a two-dimensional plot
of a complex-valued function f of a single real variable x
in the range {x0,x1}.
The plot shows the curve Abs[f] with area between the curve
and the x-axis colored by Hue[Arg[f[x]]/(2 Pi)].
The default options of Plot are changed to Axes->{True,None},
Fame->True.";

ListArgColorPlot::usage = "ListArgColorPlot[f,{x,x0,x1},opts]
plots a Abs[f], where f is a list of complex numbers. The 
points of the list Abs[f] are joined by a line.
The area between the curve and the x-axis is colored at each
point by Hue[Arg[f]/(2 Pi)].";

CombinedPlot::usage = "CombinedPlot[{f[x],g[x]},{x,x0,x1},opts]
works like ArgColorPlot with respect to f. The curve g is drawn
in front of the ArgColorPlot of f";

ListCombinedPlot::usage = "ListCombinedPlot[{list,f[x]},{x,x0,x1},opts]
works like ListArgColorPlot with respect to list.
It is assumed that list represents the discretized values of a function
defined on the interval [x0,x1]. The color list plot
is then combined with an ordinary plot of f on the same scale and
with the Ticks automatically adjusted.";

NiceTicks::usage = "NiceTicks[xmin,xmax,dx] provides a list of nice positions
for use in the Ticks or FrameTicks option in a ListPlot, where it is assumed
that the list of values ranges between xmin and xmax in steps dx.";

Saturation::usage = "Option for ArgColorPlot. Saturation->s causes the colors
in the plot to appear at saturation s. Default value is 1.";

Brightness::usage = "Option for ArgColorPlot. Brightness->b causes the colors
to be drawn at brightness b. Default value is 1.";

Begin["`Private`"]

Options[ArgColorPlot] = Join[Options[Plot],{AxesFront->True, Saturation->1,
Brightness->1}];

SetAttributes[ArgColorPlot,HoldAll]

ArgColorPlot[func_,{x_Symbol,xmin_,xmax_},opts___] :=
	Module[{comp,absfnc,argfnc,plot1,plot2,plot3,opts1,defaults,
			xvars1,xvars2,xvars3,xvars,values,hues,graph,sat,bri,style},
        {comp,sat,bri,style} =
           {Compiled,Saturation,
            Brightness,PlotStyle}/.{opts}/.Options[ArgColorPlot];
        If[comp == True,
          absfnc = Compile[{x},Abs[func]];
          argfnc = Compile[{x},
                    If[Abs[func]!=0,Mod[Arg[func]/(2 Pi),1],0]
                   ], (*else*)
          absfnc = Function[{x},N[Abs[func]]];
          argfnc = Function[{x},
                    N[If[Abs[func]!=0,Mod[Arg[func]/(2 Pi),1],0]]
                   ]
        ];
        opts1   = FilterOptions[Plot,
                  JoinOptions[DisplayFunction->Identity,opts]];
		plot1	= Plot[absfnc[x],{x,xmin,xmax}, Evaluate[opts1]];
		plot2	= Plot[ Re[func],{x,xmin,xmax}, Evaluate[opts1]];
		plot3	= Plot[ Im[func],{x,xmin,xmax}, Evaluate[opts1]];
		xvars1	= First /@ Level[plot1[[1]],{4}];Remove[plot1];
        xvars2	= First /@ Level[plot2[[1]],{4}];Remove[plot2];
        xvars3	= First /@ Level[plot3[[1]],{4}];Remove[plot3];
        xvars	= Union[xvars1,xvars2,xvars3];Remove[xvars1,xvars2,xvars3];
        values	= absfnc /@ xvars;
        hues	= Hue[argfnc[#],sat,bri]& /@ xvars;
        defaults= JoinOptions[{opts},{Axes->True, AxesFront->True}];
		Show[fillit[xvars,hues,values,style],
				Evaluate[FilterOptions[Graphics,defaults]]]
	]/;NumberQ[N[xmin]] && NumberQ[N[xmax]]


ListArgColorPlot[list_List,opts___] :=
   Module[{sat,bri,style,xvars,hues,values},
         {sat,bri,style} = {Saturation,Brightness,PlotStyle}/.{opts}/.Options[ArgColorPlot];
         xvars	= Range[Length[list]];
         hues	= Hue[Mod[Arg[#]/(2 Pi),1],sat,bri]& /@ list;
         values	= Abs[list];
         Show[fillit[xvars,hues,values,style],
             Evaluate[FilterOptions[Graphics,opts]]
         ]
   ]


CombinedPlot[{func1_,func2_},
             {x_Symbol,xmin_,xmax_},
             opts___Rule] :=
     DisplayTogether[
           ArgColorPlot[func1, {x,xmin,xmax},opts],
           Plot[func2,{x,xmin,xmax},Evaluate[FilterOptions[Plot,opts]]],
           Evaluate[FilterOptions[Graphics,opts]]
     ] /; NumberQ[N[xmin]] && NumberQ[N[xmax]]


ListCombinedPlot[{list_List,func_},
             {x_Symbol,xmin_,xmax_},
             opts___Rule] :=
   Module[{dx,x0,plotfunc},
          dx = (xmax-xmin)/(Length[list]-1);
          x0 = ind[0,xmin,dx];
          f = Function[{x},func];
          plotfunc[y_] := f[dx*(y-x0)];
          tic = genticks[xmin,xmax,dx,opts];
          DisplayTogether[
            ListArgColorPlot[list,opts],
            Plot[Evaluate[plotfunc[x]],{x,1,ind[xmax,xmin,dx]},
                Evaluate[FilterOptions[Plot,opts]]
            ],
		    FrameTicks->tic,
		    Evaluate[FilterOptions[Graphics,opts]]
		    ]
   ]/;NumberQ[N[xmin]] && NumberQ[N[xmax]]


NiceTicks[xmin_,xmax_,dx_] :=
   convertticks[LinearScale[xmin,xmax],xmin,dx] /;
   NumberQ[N[xmin]] && NumberQ[N[xmax]] && NumberQ[N[dx]]

(* auxiliary functions *)

genticks[xmin_,xmax_,dx_,opts___]:=
   Module[{tic},
          tic = FrameTicks/.{opts}/.FrameTicks->{
                               Automatic,Automatic};
          If[tic===Automatic,tic={Automatic,Automatic}];
          If[tic===None,tic={None,None}];
          If[tic[[1]]===Automatic,tic[[1]]=LinearScale[xmin,xmax]];
          If[tic[[1]]=!=None,tic[[1]] = convertticks[tic[[1]],xmin,dx]];
          tic
   ]


convertticks[list_,xmin_,dx_]:=
   Module[{res=list,a,myind},
          myind[x_] = ind[x,xmin,dx];
          Do[
             a = res[[i]];
             If[Head[a]===List,
                res[[i]]=
                If[Head[a[[2]]]===List,
                   res[[i]]=Join[{myind[a[[1]]],a[[1]]},Take[a,1-Length[a]]],
                   res[[i]]=ReplacePart[a,myind[a[[1]]],1]],
                res[[i]]={myind[a],a}
             ],
          {i,Length[res]}];
          res
   ]


ind[x_,xmin_,dx_] := N[((x-xmin)/dx)+1];


fillit[xvars_,hues_,values_,style_] :=
   Module[{nullv,xpts,valpts,lines,fills},
      nullv = Table[0,{Length[xvars]}];
      xpts = {xvars,nullv}//Transpose;
      valpts = {xvars,values}//Transpose;
      If[style === Automatic,
          lines = Line[valpts], (*else*)
          lines = Flatten[{style, Line[valpts]}]
      ];
      fills =
         {Drop[hues,-1],
             Map[Polygon,
               { Drop[xpts,-1], Drop[valpts,-1], Drop[valpts, 1], Drop[xpts, 1]
                    }//Transpose
             ]
         }//Transpose;
      Graphics[ {fills,lines} ]
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


End[]

Protect[ArgColorPlot, ListArgColorPlot, CombinedPlot, ListCombinedPlot, NiceTicks]

EndPackage[]