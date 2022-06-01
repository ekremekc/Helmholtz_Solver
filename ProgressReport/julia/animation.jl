include("fun_travwave.jl")
include("fun_Galerkin.jl")
include("fun_Helm_FD.jl")
include("fun_Helm_FE.jl")


using Plots
using LaTeXStrings


function plot_anim(number_of_elements)
TW_1, TW_2 = fun_travwave(param,1*im*pi)
GA_1, GA_2 = fun_Galerkin(param,1*im*pi, number_of_elements)
FE_1, FE_2 = fun_Helm_FE(param,1*im*pi, number_of_elements)
FD_1, FD_2 = fun_Helm_FD(param,1*im*pi, number_of_elements)

title_name = "Number of elements = " * string(number_of_elements)
title = plot(title = title_name, grid = false, showaxis = false, bottom_margin = -50Plots.px)


plot(TW_1.x, real(TW_1.P))
plot(GA_1.x, real(GA_1.P))
plot(FE_1.x[2], real(FE_1.P))
p1=plot!(FD_1.x, real(FD_1.P),legend=false,
         xlabel= L"x", ylabel= L"P_r")

plot(TW_1.x, imag(TW_1.P), label="traveling wave")
plot!(GA_1.x, imag(GA_1.P), label="Galerkin method")
plot!(FE_1.x[2], imag(FE_1.P), label="FE - Helmholtz")
p2=plot!(FD_1.x, imag(FD_1.P), label="FD - Helmholtz",
         xlabel= L"x", ylabel= L"P_{im}")

plot(TW_1.x, real(TW_1.U))
plot!(GA_1.x, real(GA_1.U))
plot!(FE_1.x[1], real(FE_1.U))
p3=plot!(FD_1.x, real(FD_1.U),legend=false,
         xlabel= L"x", ylabel= L"U_r")

plot(TW_1.x, imag(TW_1.U))
plot!(GA_1.x, imag(GA_1.U))
plot!(FE_1.x[1], imag(FE_1.U))
p4=plot!(FD_1.x, imag(FD_1.U),legend=false,
         xlabel= L"x", ylabel= L"U_{im}")



#p=plot(title,p1,p2,p3,p4,layout=(2,2))
p=plot(title,p1,p2,p3,p4,layout = @layout([A{0.05h}; [B C; D E]]))
display(p)
#savefig("my_plot.pdf")
end

n_o_el=200

anim = @animate for i ∈ 10:10:n_o_el
    plot_anim(i)
end
gif(anim, "anim_fps15.gif", fps = 6)
