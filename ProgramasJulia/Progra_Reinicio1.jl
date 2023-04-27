using LinearAlgebra, SparseArrays, Images, ImageIO, DelimitedFiles, Plots
#Método continuo------------------------------------------------------
mutable struct Particula
    radio::Real
    Posicion::Vector
    color
    Tiempo::Real
end
function Exponencial(λ)
    u = rand(1)
    m = - log.(u)./λ
    return m[1]
end
function Generacion(r, color,R,λ,n)
    t = exp(11.958211179)*n^(-0.318775)
    l1 = rand(0:360)
    y2 = [R*cos(l1*π/180), R*sin(l1*π/180)]
    τ = Exponencial(λ) * t
    return Particula(r,y2,color, τ)
end
function Pegado(Acumulado, x::Particula, y::Particula)
    Acumulado2 = copy(Acumulado)
    if minimum([norm(y.Posicion-i) for i in Acumulado]) < y.radio
             Acumulado2 = vcat(Acumulado, [y.Posicion])
             return true, Acumulado2
    end
    return false, Acumulado
#   end
end
function Radio(Acumulado)
    return maximum([norm(i) for i in Acumulado])
end

function Optimo(a::Vector,x::Particula, y::Particula)
    if norm(a + (y.Posicion - x.Posicion)/2) < x.radio
        return true
    else 
        return false
    end
end
function colision(Acumulado, x::Particula, y::Particula)
    Acumulado2 = copy(Acumulado)
    b = Acumulado[findall(a -> Optimo(a, x, y) != false, Acumulado)]
    a1 = [DSeg_Cir(i,Segmento(x.Posicion, y.Posicion)) for i in b]
    a2 = findall(x -> x != false, a1)
    if a2 != []
        a3 = minimum(a1[a2])
        a4 = Acumulado[findall(z -> DSeg_Cir(z, Segmento(x.Posicion, y.Posicion)) == a3,b)[1]]
        if a3 < y.radio
            Acumulado2 = vcat(Acumulado, [Seg_Punto(Segmento(x.Posicion, y.Posicion),a4)[2]])
            return true, Acumulado2
        end
    end
    return false, Acumulado
end
function evolucion(x::Particula)
    R = x.radio
    l1 = rand(0:360)
    y2 = [R*cos(l1*π/180), R*sin(l1*π/180)]
    return Particula(x.radio, x.Posicion + y2, x.color, x.Tiempo-1)
end
#--------------------------------------------------------------------------------------------------------------
function prueba(λ, R)
    x = Generacion(1,RGB(1.0, 0.65, 0.0), 85, λ, 2)
    Acumulado = [[0.0,0.0]]
    trac = []
    while length(Acumulado) == 1
        y = evolucion(x)
        append!(trac,[y.Posicion])
        if norm(y.Posicion) >= sqrt(2)*(R) 
            x = Generacion(x.radio, x.color,R,λ,length(Acumulado)+1)
            continue;
        end
        if y.Tiempo < 0
            break;
        end
        if norm(y.Posicion) < Radio(Acumulado) + y.radio
            if Pegado(Acumulado, x,y)[1] == true
                Acumulado = Pegado(Acumulado, x,y)[2]
                print(length(Acumulado), ",  ")
                x = Generacion(x.radio, x.color,R,λ,length(Acumulado)+1)
                continue;
            end
            if colision(Acumulado, x, y)[1] == true
                Acumulado = colision(Acumulado, x, y)[2]
                print(length(Acumulado), ",  ")
                x = Generacion(x.radio, x.color,R,λ,length(Acumulado)+1)
                continue;
            end
        end
        x = y
    end
    return trac
end
function Punto_matriz(x::Vector,n)
    L = (n+1)^2
    a = Int(floor(L/2 + x[1]))
    b = Int(floor(L/2 - x[2]))
    if a > 0 && b > 0  
      return a,b
    end
    return false
end
function Matriz_DLA1(n)
    M = spzeros((n+1)^2, (n+1)^2)
    M[Int(floor((n+1)^2/2)),Int(floor((n+1)^2/2))] = 1.0
    return M
end
function ImagenC(M)
    w, h = size(M)[1], size(M)[2]
    Imag = [RGB(0,0,0) for i in 1:w, j in 1:h]
    for i in range(1, w)
        for j in range(1,h)
            if M[i,j] != 0.0
                Imag[i,j] = RGB(1.0,1.0,1.0)
            end
        end
    end
    return Imag
end
#-------------------------------------------------------
function DLA(R,r, n,λ)
    Acumulado = [[0.0,0.0]]
    l1 = rand(0:360)
    y2 = [R*cos(l1*π/180), R*sin(l1*π/180)]
    x = Generacion(1,RGB(1.0, 0.65, 0.0), R, λ, 2)
    while length(Acumulado) < n
        y = evolucion(x)
        if norm(y.Posicion) >= sqrt(2)*(R) || y.Tiempo <= 0
            x = Generacion(x.radio, x.color,R,λ,length(Acumulado)+1)
            continue;
        end
        if norm(y.Posicion) < Radio(Acumulado) + y.radio
            if Pegado(Acumulado, x,y)[1] == true
                Acumulado = Pegado(Acumulado, x,y)[2]
                print(length(Acumulado), ",  ")
                x = Generacion(x.radio, x.color,R,λ,length(Acumulado)+1)
                continue;
            end
            if colision(Acumulado, x, y)[1] == true
                Acumulado = colision(Acumulado, x, y)[2]
                print(length(Acumulado), ",  ")
                x = Generacion(x.radio, x.color,R,λ,length(Acumulado)+1)
                continue;
            end
        end
        x = y
    end
    writedlm("Acumulado$λ.txt", Acumulado)
    return Acumulado
end
function DLA(R,r, n,λ, Acumulado)
    l1 = rand(0:360)
    y2 = [R*cos(l1*π/180), R*sin(l1*π/180)]
    x = Generacion(1,RGB(1.0, 0.65, 0.0), R, λ, 2)
    while length(Acumulado) < n
        y = evolucion(x)
        if norm(y.Posicion) >= sqrt(2)*(R) || y.Tiempo <= 0
            x = Generacion(x.radio, x.color,R,λ,length(Acumulado)+1)
            continue;
        end
        if norm(y.Posicion) < Radio(Acumulado) + y.radio
            if Pegado(Acumulado, x,y)[1] == true
                Acumulado = Pegado(Acumulado, x,y)[2]
                print(length(Acumulado), ",  ")
                x = Generacion(x.radio, x.color,R,λ,length(Acumulado)+1)
                continue;
            end
            if colision(Acumulado, x, y)[1] == true
                Acumulado = colision(Acumulado, x, y)[2]
                print(length(Acumulado), ",  ")
                x = Generacion(x.radio, x.color,R,λ,length(Acumulado)+1)
                continue;
            end
        end
        x = y
    end
    writedlm("Acumulado$λ.txt", Acumulado)
    return Acumulado
end