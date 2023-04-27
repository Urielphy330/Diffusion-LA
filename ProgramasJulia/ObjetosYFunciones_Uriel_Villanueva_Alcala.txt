mutable struct Segmento
    punto_in::Vector
    punto_fin::Vector
    longitud::Real
end

using LinearAlgebra
using Plots

function Segmento(x::Vector,y::Vector)
    Segmento(x,y,norm(x.-y))
end

mutable struct  LineaPoligonal
    vertices::Array{Vector,1}
    segmentos::Array{Segmento,1}
    longitud::Real
end

function LineaPoligonal(x)
    LineaPoligonal([i for i in x],[Segmento(x[i],x[i+1]) for i in 1:length(x)-1],sum([Segmento(x[i],x[i+1]).longitud for i in 1:length(x)-1]))
end

mutable struct Poligono
    vertices::Array{Vector,1}
    segmentos::Array{Segmento,1}
    perimetro::Real
    area::Real
end

funcion_area_del_poligono(x; y = push!(x,[x[1][1],x[1][2]])) = 1/2*abs(sum([det([y[i][1] y[i+1][1]; y[i][2] y[i+1][2]]) for i in 1:length(y)-1]))

function Poligono(x)
    Poligono([i for i in x],push!([Segmento(x[i],x[i+1]) for i in 1:length(x)-1],Segmento(x[length(x)],(x[1]))),sum([Segmento(x[i],x[i+1]).longitud for i in 1:length(x)-1])+Segmento(x[length(x)],(x[1])).longitud,funcion_area_del_poligono([[x[i][1],x[i][2]] for i in 1:length(x)]))
end

mutable struct SuperficiePoliedrica
    poligonos::Array{Poligono,1}
    superficie::Real
end

function SuperficiePoliedrica(x)
    SuperficiePoliedrica([i for i in x],sum([i.Area for i in x]))
end

mutable struct Poliedro
    caras::Array{Poligono,1}
    superficie::Real
    volumen::Real
end

import Base.+

function +(x::Segmento,y::Vector)
    if length(x.punto_fin) == length(y)
        return LineaPoligonal([x.punto_in,x.punto_fin,y])
    else
        return(error("Error de dimensión"))
    end
end

function +(y::Vector, x::Segmento)
    return +(x,y)
end

function +(x::Segmento, y::Segmento)
    if x.punto_fin == y.punto_in 
        return LineaPoligonal([x.punto_in,x.punto_fin,y.punto_fin])
    else 
        return(error("Error de posición"))
    end
end

function +(x::LineaPoligonal, y::Vector)
    if length(x.vertices[1]) == length(Tuple(y))
        return LineaPoligonal(push!([i for i in x.vertices],y))
    else
        return error("Error de dimensión")
    end
end

function +(y::Vector, x::LineaPoligonal)
    return +(x,y)
end

function +(x::LineaPoligonal, y::Segmento)
    if x.vertices[end] == y.punto_in
        return LineaPoligonal(push!([i for i in x.vertices],y.punto_fin))
    else
        return error("Error de posición")
    end
end

function +(y::Segmento,x::LineaPoligonal)
    return +(x,y)
end

function +(x::LineaPoligonal,y::LineaPoligonal)
    if x.vertices[end] == y.vertices[1]
        return LineaPoligonal(append!([i for i in x.vertices],[y.vertices[i] for i in 2:length(x.vertices)]))
    else
        return error("Error de posición")
    end
end

function +(x::Poligono, y::Poligono)
    if sum([i==j for i in x.segmentos, j in y.segmentos]) == 1
        return SuperficiePoliedrica([x,y], sum(x.Area+y.Area))
    else
        return error("Error")
    end
end

function +(x::SuperficiePoliedrica, y::Poligono)
    for k in x.poligonos
        if sum([i == j for i in k.segmentos, j in y.segmentos]) >= 1
            return SuperficiePoliedrica(push!([i for i in x.poligonos],y))
        end
    end
    return error("Error")
end

function +(y::Poligono, x::SuperficiePoliedrica)
    return +(x,y)
end

function +(x::SuperficiePoliedrica, y::SuperficiePoliedrica)
    for k in x.poligonos
        for l in y.poligonos
            if sum([i == j for i in k.segmentos, j in l.segmentos]) >= 1
                return SuperficiePoliedrica(append!([i for i in x.poligonos],[j for j in y.poligonos]))
            end
        end
    end
    return error("Error")
end
function pertenencia(x::Segmento, y::Vector) #version correcta
    if x.punto_fin == y || x.punto_in == y
        return true
    end
    if x.punto_fin[2] == x.punto_in[2] 
        return 0 < (y[1]-x.punto_in[1])/(x.punto_fin[1]-x.punto_in[1]) < 1
    end
    t = (y[2] - x.punto_in[2])/(x.punto_fin[2]-x.punto_in[2])
    return (x.punto_fin[2]-x.punto_in[2])/(x.punto_fin[1]-x.punto_in[1]) ≈ (y[2]-x.punto_in[2])/(y[1]-x.punto_in[1]) && 0<t<1  
end
struct Tiro_parabolico
    v::Vector
    p::Vector
end
function tiro_parabolico(x::Tiro_parabolico; n = 1000)
    t_0 = x.v[2]/9.81 + sqrt(x.v[2]^2+19.62*x.p[2])/9.81
    return [(x.p[1]+x.v[1]*t,x.p[2]+x.v[2]*t-1/2*9.81*t^2) for t in 0:t_0/n:t_0]
end
struct Recta
    f::Vector
    m::Real
end
function recta(x::Recta)
    if x.m == Inf || x.m == -Inf
        [(x.f[1],i) for i in -2:0.1:10]
    else
        [(i,x.m*(i-x.f[1])+x.f[2]) for i in -2:0.1:15]
    end
end
function inter_tiro_rec(x::Tiro_parabolico,y::Recta) #Version correcta
    if y.m == Inf || y.m == -Inf
       if (y.f[1]-x.p[1])/x.v[1] < (-y.v[2]-sqrt((x.v[2])^2-4(-4.05)(x.p[2])))/-9.81
            return true, [(y.f[1]-x.p[1])/x.v[1],x.p[2]]
       else
            return false
       end
    end
    if ((x.v[2]-y.m*x.v[2])^2-4(-4.905)(x.p[2]-y.f[2]-y.m*(x.p[1]-y.f[1]))) < 0
        return false
    else
    t_1 = -(x.v[2]-y.m*x.v[1])/-9.81 + sqrt((x.v[2]-y.m*x.v[1])^2-4(-4.905)(x.p[2]-y.f[2]-y.m*(x.p[1]-y.f[1])))/-9.81
    t_2 = -(x.v[2]-y.m*x.v[1])/-9.81 - sqrt((x.v[2]-y.m*x.v[1])^2-4(-4.905)(x.p[2]-y.f[2]-y.m*(x.p[1]-y.f[1])))/-9.81
    x_1 = x.p[1]+x.v[1]*t_1
    x_2 = x.p[1]+x.v[1]*t_2
    y_1 = y.m*(x_1-y.f[1])+y.f[2]
    y_2 = y.m*(x_2-y.f[1])+y.f[2]
    return [x_1,y_1,t_1],[x_2,y_2,t_2]
    end
end
function inter_seg_tiro(x::Tiro_parabolico, y::Segmento) #version correcta
    m = (y.punto_fin[2]-y.punto_in[2])/(y.punto_fin[1]-y.punto_in[1])
    a = Recta(y.punto_in, m)
    b = inter_tiro_rec(x,a)
    if b == false
       return false
    end
    if pertenencia(y,[b[1][1],b[1][2]]) == true && pertenencia(y,[b[2][1],b[2][2]])==true
      return [b[1][1],b[1][2],b[1][3]], [b[2][1],b[2][2],b[2][3]]
    else
       if pertenencia(y,[b[1][1],b[1][2]]) == true
           return [b[1][1],b[1][2],b[1][3]]
       end
       if pertenencia(y,[b[2][1],b[2][2]]) == true
           return [b[2][1],b[2][2],b[2][3]]
        end
        return false
    end
end
function inter_tiro_linea(x::Tiro_parabolico,y::LineaPoligonal) #version correcta
    h = [inter_seg_tiro(x,i) for i in y.segmentos]
    if sum([i != false for i in h]) == 0 
        return false
    else
        a = findall(x -> x != false,h)
        if length(a) == 0
            return false
        end
        b = sort(h[a], by = x-> x[3])
        return [b[1][1],b[1][2],b[1][3]]
    end
end
                    
                    
import Plots.plot
import Plots.plot!
function plot(x::Segmento; karg...)
    plot([x.punto_in[1],x.punto_fin[1]],[x.punto_in[2],x.punto_fin[2]];karg...)
    scatter!([Tuple(x.punto_in), Tuple(x.punto_fin)], color = :blue, key = false, markersize = 0.8)
end
function plot!(x::Segmento; karg...)
    plot!([x.punto_in[1],x.punto_fin[1]],[x.punto_in[2],x.punto_fin[2]];karg...)
    scatter!([Tuple(x.punto_in), Tuple(x.punto_fin)], color = :blue, key = false, markersize = 0.8)
end
function plot(x::LineaPoligonal; karg...)
    plot(x.segmentos[1]; karg...)
    for i in 2:length(x.segmentos)
        plot!(x.segmentos[i]; karg...)
    end
    plot!()
end
function plot!(x::LineaPoligonal; karg...)
    for i in x.segmentos
        plot!(i; karg...)
    end
    plot!()
end
function plot(x::Poligono; karg...)
    plot(x.segmentos[1]; karg...)
    for i in 2:length(x.segmentos)
        plot!(x.segmentos[i]; karg...)
    end
    plot!()
end
function plot!(x::Poligono; karg...)
    for i in x.segmentos
        plot!(i; karg...)
    end
    plot!()
end
function plotcirculo!(p::Vector,r::Real; karg...)
    θ = 0:0.05:(2*π + 0.05)
    x = r*(cos.(θ) .+ p[1])
    y = r*(sin.(θ) .+ p[2])
    plot!(x,y,aspect_ratio = :equal; karg...)
end
  
function Seg_Punto(a::Segmento, x::Vector)
    m = (a.punto_fin[2] - a.punto_in[2])/(a.punto_fin[1] - a.punto_in[1])
    if m == 0
        return pertenencia(a, [x[1],a.punto_in[2]]), [x[1],a.punto_in[2]]
    elseif m == Inf
        return pertenencia(a, [a.punto_in[1], x[2]]), [a.punto_in[1], x[2]] 
    else
    x1 = (m^(-1)*x[1] + x[2] - a.punto_in[2]+m*a.punto_in[1])/(m + m^(-1)) 
    y1 = x[2] - m^(-1)*(x1 - x[1])
    return pertenencia(a, [x1,y1]), [x1, y1]
    end
end
function DSeg_Cir(c::Vector, a::Segmento)
    if Seg_Punto(a,c)[1] != false
        return norm(Seg_Punto(a,c)[2]-c)
    else
        return false
    end
end
function Seg_Punto(a::Segmento, x::Vector)
    m = (a.punto_fin[2] - a.punto_in[2])/(a.punto_fin[1] - a.punto_in[1])
    if m == 0
        return pertenencia(a, [x[1],a.punto_in[2]]), [x[1],a.punto_in[2]]
    elseif m == Inf
        return pertenencia(a, [a.punto_in[1], x[2]]), [a.punto_in[1], x[2]] 
    else
    x1 = (m^(-1)*x[1] + x[2] - a.punto_in[2]+m*a.punto_in[1])/(m + m^(-1)) 
    y1 = x[2] - m^(-1)*(x1 - x[1])
    return pertenencia(a, [x1,y1]), [x1, y1]
    end
end
function DSeg_Cir(c::Vector, a::Segmento)
    if Seg_Punto(a,c)[1] != false
        return norm(Seg_Punto(a,c)[2]-c)
    else
        return false
    end
end
function diametro(x::Poligono)
    return(maximum([norm(i-j) for i in x.vertices,j in x.vertices]))
end
function interseccion_segmentos(x::Segmento, y::Segmento)
    if x == y 
        return Inf
    end
    m1 = (x.punto_fin[2]-x.punto_in[2])/(x.punto_fin[1]-x.punto_in[1])
    m2 = (y.punto_fin[2]-y.punto_in[2])/(y.punto_fin[1]-y.punto_in[1])
    if  m1 == m2 
        return false
    end
    A = [y.punto_fin[1]-y.punto_in[1] -(x.punto_fin[1]-x.punto_in[1]);
         y.punto_fin[2]-y.punto_in[2] -(x.punto_fin[2]-x.punto_in[2])]
    b = [x.punto_in[1]-y.punto_in[1]; x.punto_in[2]-y.punto_in[2]]
    t = A\b
    if t[1] < 0 || t[2] < 0 || t[1] > 1  || t[2] > 1
        return false
    else
        return t[1]*(y.punto_fin-y.punto_in)+y.punto_in
    end     
end
function dentro_fuera(x::Poligono,y::Vector)
    m = rand()
    a = Segmento(y,[sqrt(diametro(x)^2/(1 + m^2))+y[1], m * sqrt(diametro(x)^2/(1+m^2))+y[2]])
    b = sum([interseccion_segmentos(a,i) != false for i in x.segmentos])
    if b % 2 == 0
        return false
    else
        return true
    end
end