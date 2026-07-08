using Revise, ADerrors, FitRoutines, PyPlot

xdata = collect(0:10)
true_model(x) = 1.0 + 3exp(-x/5)

ydata = let
    y = true_model.(xdata)
    [uwreal(randn(1000).+_y, "test") for _y in y]
end

Q = uwreal([10.0,10.0], "prior")
uwerr(Q)
model(x,p,q) = @. p[1] + p[2]*exp(-x/q[1])

model(x,p)= p[1] + p[2]*exp(-x/p[3])

fit1 = FitRoutines.fit_routine_prior(model,xdata,ydata,[Q],2,corr=false)
fit2 = FitRoutines.fit_routine_prior(model,xdata,ydata,[(Q.mean,Q.err)],
                                     2,corr=false)
0==0
