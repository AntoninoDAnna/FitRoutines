using Revise, ADerrors, FitRoutines, PyPlot

xdata = [collect(0:10) for _ in 1:2]
true_model1(x) = 1.0 + 3exp(-x/5)
true_model2(x) = 1.0 - 0.5*exp(-x/5)

ydata1 = let
    y = true_model1.(xdata[1])
    [uwreal(randn(1000).+_y, "test") for _y in y]
end

ydata2 = let
    y = true_model2.(xdata[2])
    [uwreal(randn(1000).+_y, "test") for _y in y]
end
ydata = [ydata1,ydata2]

Q = uwreal([10.0,10.0], "prior")
uwerr(Q)

model1(x,p,q) = @. p[1] + p[2]*exp(-x/q[1])
model2(x,p,q) = @. p[1] + p[3]*exp(-x/q[1])

model1(x,p)=  p[1] + p[2]*exp(-x/p[4])
model2(x,p)=  p[1] + p[3]*exp(-x/p[4])

model = [model1,model2]

fit1 = FitRoutines.fit_routine_prior(model,xdata,ydata,[Q],3,corr=false)
fit2 = FitRoutines.fit_routine_prior(model,xdata,ydata,[(Q.mean,Q.err)],
                                     3,corr=false)
0==0


let
    fig,ax = subplots(1)
    for (x,y,c,m,mod) in zip(xdata,ydata,["k","r"],["x","s"],model)
        uwerr.(y)
        ax.errorbar(x,value.(y),ADerrors.err.(y), color=c, marker=m,
                    fillstyle="none",ls="")
        _y =[mod(_x,fit1.par) for _x in x]
        uwerr.(_y)
        v,e = value.(_y),ADerrors.err.(_y)
        ax.fill_between(x,v.-e,v.+e, color="g")
    end
end
