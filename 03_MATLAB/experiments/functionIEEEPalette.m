function S = functionIEEEPalette()
%FUNCTIONIEEEPALETTE  Shared figure vocabulary for the Paper 2 figure set.
%
%Returns a struct of colours and ready-made series specifications so that every
%figure in the manuscript uses the same encoding. The design rule is that NO
%series is distinguished by colour alone: each carries a colour, a line style
%and a marker together, so the figure survives a greyscale print and a
%colour-blind reader. This is the IEEE expectation for figures that will be
%printed as well as read on screen.
%
%  S.ink / S.brick / S.teal / S.ochre / S.grey   RGB triples
%  S.s1 ... S.s4      series specs, each a struct with .c .ls .mk
%  S.fs*              the type sizes used throughout
%
%Sizes assume the figure is exported at its true printed width (8.8 cm for an
%IEEE single column), so the point sizes set here are the point sizes the
%reader sees. Do not scale the figure after export.
%
%Version 1.1 (2026-09-19). License: GPLv2.

S.ink   = [0.102 0.102 0.102];      % #1A1A1A
S.brick = [0.706 0.263 0.169];      % #B4432B
S.teal  = [0.122 0.435 0.471];      % #1F6F78
S.ochre = [0.753 0.541 0.118];      % #C08A1E
S.grey  = [0.550 0.550 0.550];
S.faint = [0.870 0.870 0.870];

%series 1..4, in the order they should be consumed
S.s1 = struct('c',S.ink  ,'ls','-'  ,'mk','o');
S.s2 = struct('c',S.brick,'ls','--' ,'mk','s');
S.s3 = struct('c',S.teal ,'ls','-.' ,'mk','^');
S.s4 = struct('c',S.ochre,'ls',':'  ,'mk','d');

S.fnt   = 'Times New Roman';        % matches IEEEtran body type
S.fsTick = 8;
S.fsLab  = 9;
S.fsLeg  = 8.0;   % IEEE floor for legend type at final print size; the widest
                  % legend strings were shortened to fit rather than shrinking type
S.fsAnn  = 7;

S.lw     = 1.1;                     % data line width
S.lwAx   = 0.7;                     % axis and box line width
S.ms     = 4.5;                     % marker size
S.cap    = 2.5;                     % error-bar cap size

S.wCol   = 8.8;                     % IEEE single-column width [cm]
end
