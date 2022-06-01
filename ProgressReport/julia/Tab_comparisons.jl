include("fun_travwave.jl")
include("fun_Galerkin.jl")
include("fun_Helm_FD.jl")
include("fun_Helm_FE.jl")


using Plots
using LaTeXStrings

number_of_elements=400
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
savefig("my_plot.pdf")


## GENERATE TABLE AND WRITE INTO EXCEL FILE
using CSV, DataFrames

table = ["N" "fun_travwave" "fun_Helm_FD" "fun_Helm_FE" "fun_Galerkin";
         1 TW_1.s "-" "-" "-";
         ]

N_array = [10,40,100,400]

for N in N_array
        gal_temp,gal_garb =fun_Galerkin(param,1*im*pi, N)
        FD_temp, FD_garb = fun_Helm_FD(param,1*im*pi,N)
        FE_temp, FE_garb = fun_Helm_FE(param,1*im*pi,N)
        new_row = [N "-" gal_temp.s FD_temp.s FE_temp.s]
        global table = [table;new_row]
end

CSV.write("eigenvalue_table.csv", DataFrame(table),header=false)

## GENERATE SENSTIVITIES TABLE AND WRITE INTO EXCEL FILE


gal_temp,gal_garb =fun_Galerkin(param,1*im*pi, number_of_elements)
TW_1, TW_2 = fun_travwave(param,1*im*pi)
FD_temp, FD_garb = fun_Helm_FD(param,1*im*pi,number_of_elements)
FE_temp, FE_garb = fun_Helm_FE(param,1*im*pi,number_of_elements)

table_sens = [" " "∂s/∂n" "∂s/∂τ" "∂s/∂x_m" "∂s/∂x_f";
              "fun_travwave.jl" TW_2.n TW_2.tau TW_2.x_m TW_2.x_f;
              "fun_Helm_FD.jl" FD_garb.n FD_garb.tau FD_garb.x_m FD_garb.x_f;
              "fun_Helm_FE.jl" FE_garb.n FE_garb.tau FE_garb.x_m FE_garb.x_f;
              "fun_Galerkin.jl" gal_garb.n gal_garb.tau gal_garb.x_m gal_garb.x_f]



CSV.write("senstivities_table.csv", DataFrame(table_sens),header=false)
