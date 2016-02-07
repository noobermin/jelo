<TeXmacs|1.99.2>

<style|generic>

<\body>
  For a constant magnetic field,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|d<around*|(|\<gamma\><wide|\<beta\>|\<vect\>>|)>|
    d t>>|<cell|=>|<cell|<frac|q B|m><wide|\<beta\>|\<vect\>>\<times\><wide|z|^>.>>>>
  </eqnarray*>

  The expression on the left-hand side has units of inverse time. Thus, we
  have that

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|[|<frac|q
    B|m>|]>>|<cell|=>|<cell|T<rsup|-1>,>>>>
  </eqnarray*>

  and <math|q B/m> is referred to as the magnetic gyrofequency
  <math|\<omega\><rsub|g>>, or cylcotron frequency, in some contexts.
  Therefore,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|d<around*|(|\<gamma\><wide|\<beta\>|\<vect\>>|)>|d<around*|(|\<omega\><rsub|g>t|)>>>|<cell|=>|<cell|\<beta\><rsub|y><wide|x|^>-\<beta\><rsub|x><wide|y|^>.>>>>
  </eqnarray*>

  In terms of the differential proper time <math|d u =d t/\<gamma\>>,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|d<around*|(|\<gamma\><wide|\<beta\>|\<vect\>>|)>|d
    u>>|<cell|=>|<cell|<around*|[|\<gamma\>\<beta\><rsub|y><wide|x|^>-\<gamma\>\<beta\><rsub|x><wide|y|^>|]>\<omega\><rsub|g>>>|<row|<cell|\<Rightarrow\><frac|d<around*|(|\<gamma\><wide|\<beta\>|\<vect\>>|)>|d
    u>\<cdot\>\<gamma\><wide|\<beta\>|\<vect\>>>|<cell|=>|<cell|<frac|d\<gamma\>|d
    u>>>|<row|<cell|>|<cell|=>|<cell|0>>>>
  </eqnarray*>

  which implies that the gamma factor is a constant of motion. Thus,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|d<around*|(|\<gamma\><wide|\<beta\>|\<vect\>>|)>|d
    u>>|<cell|=>|<cell|\<gamma\><frac|d<wide|\<beta\>|\<vect\>>|d
    u>>>|<row|<cell|>|<cell|=>|<cell|<around*|[|\<gamma\>\<beta\><rsub|y><wide|x|^>-\<gamma\>\<beta\><rsub|x><wide|y|^>|]>\<omega\><rsub|g>>>|<row|<cell|\<Rightarrow\><frac|d<wide|\<beta\>|\<vect\>>|d
    u>>|<cell|=>|<cell|\<omega\><rsub|g><around*|[|\<beta\><rsub|y><wide|x|^>-\<beta\><rsub|x><wide|y|^>|]>>>|<row|<cell|\<Rightarrow\><frac|d
    \<beta\><rsub|x>|d u>>|<cell|=>|<cell|\<omega\><rsub|g>\<beta\><rsub|y>,>>|<row|<cell|<frac|d
    \<beta\><rsub|y>|d u>>|<cell|=>|<cell|-\<omega\><rsub|g>\<beta\><rsub|x>>>|<row|<cell|\<Rightarrow\><frac|d
    <rsup|2>\<beta\><rsub|y>|d u<rsup|2>>>|<cell|=>|<cell|-\<omega\><rsub|g><rsup|2>\<beta\><rsub|y>>>|<row|<cell|\<Rightarrow\>\<beta\><rsub|y><around*|(|u|)>>|<cell|=>|<cell|\<beta\><rsub|c>cos<around*|(|\<omega\><rsub|g>u|)>+\<beta\><rsub|s>sin<around*|(|\<omega\><rsub|g>u|)>,>>|<row|<cell|\<beta\><rsub|x><around*|(|u|)>>|<cell|=>|<cell|\<beta\><rsub|c>sin<around*|(|\<omega\><rsub|g>u|)>-\<beta\><rsub|s>cos<around*|(|\<omega\><rsub|g>u|)>.<eq-number><label|sol>>>>>
  </eqnarray*>

  with constants fixed by boundary conditions. Specifically, we find that

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|int><rsub|0><rsup|t>d
    t<rprime|'>>|<cell|=>|<cell|<big|int><rsub|u<rsub|0>><rsup|u>\<gamma\>d
    u<rprime|'>>>|<row|<cell|t>|<cell|=>|<cell|\<gamma\>\<Delta\>u.>>>>
  </eqnarray*>

  Defining <math|u<rsub|0>=u<around*|(|t=0|)>>, we have that
  <math|t=\<gamma\><around*|(|u-u<rsub|0>|)>>. Thus, absorbing
  <math|u<rsub|0>> into the constants in (<reference|sol>), we obtain

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|\<beta\>|\<vect\>><around*|(|t|)>>|<cell|=>|<cell|<matrix|<tformat|<table|<row|<cell|\<beta\><rsub|c>sin<around*|(|\<omega\><rsub|g>t|)>-\<beta\><rsub|s>cos<around*|(|\<omega\><rsub|g>t|)>>>|<row|<cell|\<beta\><rsub|c>cos<around*|(|\<omega\><rsub|g>t|)>+\<beta\><rsub|s>sin<around*|(|\<omega\><rsub|g>t|)>>>|<row|<cell|\<beta\><rsub|z0>>>>>>>>>>
  </eqnarray*>

  where we defined the ``relativistic gyrofrequency''
  <math|\<omega\><rsub|g><rprime|'>=\<omega\><rsub|g>/\<gamma\>=<frac|q
  B|\<gamma\>m>.> To specify <math|<wide|\<beta\>|\<vect\>><around*|(|0|)>=<wide|\<beta\>|\<vect\>><rsub|0>>
  we can fix that

  <\eqnarray*>
    <tformat|<table|<row|<cell|<matrix|<tformat|<table|<row|<cell|\<beta\><rsub|x0>>>|<row|<cell|\<beta\><rsub|y0>>>|<row|<cell|\<beta\><rsub|z0>>>>>>>|<cell|=>|<cell|<matrix|<tformat|<table|<row|<cell|-\<beta\><rsub|s>>>|<row|<cell|\<beta\><rsub|c>>>|<row|<cell|\<beta\><rsub|z0>>>>>>>>|<row|<cell|\<Rightarrow\><wide|\<beta\>|\<vect\>><around*|(|t|)>>|<cell|=>|<cell|<matrix|<tformat|<table|<row|<cell|\<beta\><rsub|y0>sin<around*|(|\<omega\><rsub|g>t|)>+\<beta\><rsub|x0>cos<around*|(|\<omega\><rsub|g>t|)>>>|<row|<cell|\<beta\><rsub|y0>cos<around*|(|\<omega\><rsub|g>t|)>-\<beta\><rsub|x0>sin<around*|(|\<omega\><rsub|g>t|)>>>|<row|<cell|\<beta\><rsub|z0>>>>>>>>>>
  </eqnarray*>

  From this, one can see that <math|\<beta\><rsup|2><around*|(|t|)>=<left|\|><wide|\<beta\>|\<vect\>><rsub|0><left|\|><rsup|2>>,
  so that <math|\<gamma\>> is a constant of motion in a particular frame. For
  the sake of completeness, Thus, the full motion is

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|x|\<vect\>><around*|(|t|)>>|<cell|=>|<cell|<wide|x|\<vect\>><rsub|0>+<frac|1|\<omega\><rsub|g>><matrix|<tformat|<table|<row|<cell|\<beta\><rsub|y<rsub|0>><around*|(|1-cos<around*|(|\<omega\><rsub|g>t|)>|)>+\<beta\><rsub|x0>sin<around*|(|\<omega\><rsub|g>t|)>>>|<row|<cell|\<beta\><rsub|y0>sin<around*|(|\<omega\><rsub|g>t|)>+\<beta\><rsub|x0><around*|(|cos<around*|(|\<omega\><rsub|g>t|)>-1|)>>>|<row|<cell|\<beta\><rsub|z0>\<omega\><rsub|g>t>>>>>.>>>>
  </eqnarray*>
</body>

<\references>
  <\collection>
    <associate|sol|<tuple|1|?|.TeXmacs/texts/scratch/no_name_5.tm>>
  </collection>
</references>