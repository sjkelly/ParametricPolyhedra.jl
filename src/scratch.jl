using SymPy

a = Sym(:a)
b = Sym(:b)
c = Sym(:c)
aa = Sym(:aa)
ba = Sym(:ba)
ca = Sym(:ca)

solve([a*cos(ba) + b*cos(aa) - c,
       b*sin(aa) - a*sin(ba),
       aa + ba + ca - pi],[aa,ba,ca,a])

@profile solve([a*cos(ba) + b*cos(aa) - c,
       b*sin(aa) - a*sin(ba),
       aa + ba + ca - pi],[aa,ba,ca,a])

using ProfileView
ProfileView.view()

