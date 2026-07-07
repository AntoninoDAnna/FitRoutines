using Revise, ADerrors, FitRoutines, PyPlot

xdata = collect(1:10)
true_model(x) = 1.0 + 3exp(-x/5)
ydata = let
    y = true_model.(xdata)
    [uwreal(randn(1000).+_y, "test") for _y in y]
end

Q = uwreal([5.,0.1], "prior")

model(x,p,q) = @. p[1] + p[2]*exp(-x/q[1])
model(x,p)= p[1] + p[2]*exp(-x/p[3])
fit = FitRoutines.fit_routine_prior(model,xdata,ydata,[Q],2,corr=true)

let
    fig,ax =subplots(1)
    y,dy = value.(ydata), ADerrors.err.(ydata)
    ax.errorbar(xdata,y,dy,color="k",marker="s",fillstyle="none",ls="")
    r = collect(1:0.1:10)
    ax.plot(r,true_model.(r),color="r",ls="-")
    y = [model(x,fit.par) for x in r]
    uwerr.(y)
    y,dy = value.(y), ADerrors.err.(y)
    ax.fill_between(r,y.-dy, y.+dy,color="g",alpha=0.5)
end
