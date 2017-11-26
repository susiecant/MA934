function keyvalprint(list::Nullable{LList})
    #=This function traverses inputted list and prints out the key-value pairs
    stored in it=#
    if isnull(list)
    return           #terminates at end of list

    else
    println(get(list).data)
    keyvalprint(get(list).next)        #recursive line
    end
end
function searchllist(list::Nullable{LList},k::Int64)
#=This function searches an LList for the key k and returns the corresponding
KVPair if it is present and a Nullable{KVPair} otherwise=#
        if isnull(list)
        return Nullable{KVPair}

        elseif (get(list).data.key==k)
        return get(list).data

        else
        searchllist(get(list).next,k)   #recursive line
        end
end
function intervalmembership(list::Nullable{LList},x::Float64)
#=This function takes the LList containing the list of partial sums and a random
Float64 in the range [0,xn] as inputs and returns the KVPair corresponding
to the interval in which x lies. =#
        if isnull(list)
            return

        elseif (x<get(list).data.value)
            return get(list).data

        else
        intervalmembership(get(list).next,x)     #recursive line
        end
end
function intervalmembership(FT::Nullable{FTree},x::Float64)
#=This function takes the array of KVPairs containing the interval lengths as
input, recursively constructs the tree and returns the FTree containing the
correct key-value pairs=#
        if (isnull(get(FT).left) && isnull(get(FT).right))
            return get(FT).data

        elseif (x<get(get(FT).left).data.value)
            intervalmembership(get(FT).left,x)

        else
            intervalmembership(get(FT).right,x-get(get(FT).left).data.value)
            #recursive line
        end
    end

function Gillespie_list(N)
#=This function outputs X a range of positions and P the particle density at
each of these positions after running a Gillespie algorithm with random
exponential raters.=#
L=10.0
Nx = 201
dx = 2.0*L/(Nx-1)
X = dx.*(-(Nx-1)/2:(Nx-1)/2)
Y =zeros(Int64,N)
t=0.0
T=1.0

D=randexp(N)
D=vcat(D,D)          #same rate to move left or right
r = (D./2.0)/(dx*dx)
totalRate = sum(r)
dt = 1.0/totalRate
cum_r=cumsum(r)
#cumulative sum of rates to give coorect structure for LList search
values = Array{KVPair}(2*N)
for i in 1:2*N
    values[i] = KVPair(i,cum_r[i])
end
myList=Nullable{LList}()  #build linked list
myList=buildLList(values)

# This is the main loop
while t < T
    rand_k=rand()*totalRate #random number between [0,R]
    k=intervalmembership(myList,rand_k).key
    #search linked list for the interval containing random number
    if k<=N
        hop = 1
        particleId = k

    else
        hop = -1
        particleId=k-N

    end

    Y[particleId]+=hop
    t+=dt
end


# Calculate the estimated density of particles
P =zeros(Float64,length(X))
for i in 1:length(Y)
    P[Y[i]+Int64((Nx-1)/2)+1]+=1/(N * dx)
end
    return X,P
end

function Gillespie_tree(N)
    #=This function outputs X a range of positions and P the particle density at
    each of these positions after running a Gillespie algorithm with random
    exponential raters.=#
L=10.0
Nx = 201
dx = 2.0*L/(Nx-1)
X = dx.*(-(Nx-1)/2:(Nx-1)/2)
Y =zeros(Int64,N)
t=0.0


D=randexp(N)
D=vcat(D,D)    #same rate to move left and right
r = (D./2.0)/(dx*dx)
totalRate = sum(r)
dt = 1.0/totalRate
T=1.0
values = Array{KVPair}(2*N)
for i in 1:2*N
    values[i] = KVPair(i,r[i])
end
myTree = Nullable{FTree}(FTree(KVPair(0,0.0))) #Build tree
myTree=buildFTree(myTree, values);

# This is the main loop
while t < T
    rand_k=rand()*totalRate #random number between [0,R]
    k=intervalmembership(myTree,rand_k).key
    #find which interval random number lies in
    if k<=N
        hop = 1
        particleId = k

    else
        hop = -1
        particleId=k-N

    end

    Y[particleId]+=hop
    t+=dt
end


# Calculate the estimated density of particles
P =zeros(Float64,length(X))
for i in 1:length(Y)
    P[Y[i]+Int64((Nx-1)/2)+1]+=1/(N * dx)
end
return X,P
end
function diffusionpdf(x, D, t)
    return (1.0/sqrt(2.0*D*t))*exp(-abs(x)*sqrt(2/(D*t)))
end
