function finitediff3(func,x,h)
#Gives the 3 point stencil finite difference formula for the derivative of func
    return (4*func(x+h)-3*func(x)-func(x+2*h))/(2*h)
end
function sinexp(x)
#Evaluates sin(exp(x)) at x
    return sin(exp(x))
end
function Dsinexp(x)
#Evaluates the derivative of sin(exp(x))
    return exp(x)*cos(exp(x))
end
function bracketbisect(func,a,b,epstol)
    error=zeros(100)
    i=1
    while abs(b-a)>epstol
        if func(a)*func((b+a)/2)<0
            a=a
            b=(b+a)/2
            error[i]=abs(((b+a)/2)-log(pi))
            i=i+1
        else
            a=(b+a)/2
            b=b
            error[i]=abs(((b+a)/2)-log(pi))
            i=i+1
        end
        end
    return (b+a)/2,error[1:i-1],i-1
end
function NewtonRaph(func,Dfunc,x,eps)
    error=zeros(100)
    i=1
    while abs(func(x))>eps
        x_new=x-(func(x)/Dfunc(x))
        x=x_new
        error[i]=abs(x-log(pi))
        i=i+1
    end
    return x,error[1:i-1],i
end
function Goldensec(func,a,b,epstol)
    phi=(1+sqrt(5))/2
    c = b-((b-a)/phi)
    error=zeros(100)
    i=1
    while abs(b-a)>epstol
    x = a+((b-a)/phi)
    c = b-((b-a)/phi)

    if func(x)<func(c)
            error[i]=abs(c-log(1.5*pi))
            i=i+1
            a=c
            c=x
            b=b

        else
            error[i]=abs(c-log(1.5*pi))
            i=i+1
            a=a
            c=c
            b=x

        end
    end
    return c,error[1:i-1],i-1
end
