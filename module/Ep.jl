# Insert absolute path of "Initialize" module
push!(LOAD_PATH, "/home/sl18/julia_scripts/julia_master/Speedcat/module/")

# Define module
module Ep

# Make it visible outside
export solve_ep

# Import modules
using Initialize
using Revise
using Statistics
using Sundials

# Gattoni rat heart EP model (2016)
function fun(du, u, p, t)
	R, T, F, Cm, stim_period, stim_duration, stim_amplitude, Vmyo, Vmyo_uL, VSR_uL, g_Na, Na_o, g_to, a, b, K_o , g_ss , g_K1, g_f, f_Na, g_B_Na, g_B_K, i_NaK_max, gamma1, KdNaio, KdNais, KdNaes, Delta, alpha, Ca_o, g_D, J_R, J_L, N, V_L, del_VL, phi_L, t_L, tau_L, t_R, tau_R, phi_R, theta_R, K_RyR, K_L, a1, b1, c, d, K_mNa, K_mCa, eta, k_sat, g_NCX, g_SERCA, K_SERCA , g_pCa, K_mpCa, g_CaB, g_SRl, k_m_TRPN, k_p_TRPN, B_TRPN, k_CMDN, B_CMDN, EGTA_tot, KmEGTA, tau_s_ss, f_K, alpha_m, beta_m = p
	z_1, z_2, z_3, r, s, s_slow, y, Ca_SR, Ca_i, K_i, Na_i, TRPN, V, h, j, m, r_ss, s_ss = u

	i_Stim = (t-fld(t, stim_period)*stim_period >= 0.0 && t-fld(t, stim_period)*stim_period <= stim_duration) ? stim_amplitude : 0.0

	FVRT = F*V/(R*T)
	FVRT_Ca = 2.0*FVRT

	E_Na = R*T/F*log(Na_o/Na_i)
	i_Na = g_Na*m^3*h*j*(V-E_Na)

	m_inf = (1.0+exp((V+45.0)/(-6.5)))^(-1)
	tau_m = 1.36/(((0.32*(V+47.13))/(1.0-exp(-0.1*(V+47.13))))+0.08*exp(-V/11.0))

	h_inf = (1.0+exp((V+76.1)/6.07))^(-1)
	tau_h = (V >= -40.0) ? 0.4537*(1.0+exp(-(V+10.66)/11.1)) : 3.49/(0.135*exp(-(V+80.0)/6.8)+3.56*exp(0.079*V)+310000.0*exp(0.35*V))

	j_inf = (1.0+exp((V+76.1)/6.07))^(-1)
	tau_j = (V >= -40.0) ? 11.63*(1.0+exp(-0.1*(V+32.0)))*exp(0.0000002535*V) : 3.49/((V+37.78)/(1.0+exp(0.311*(V+79.23)))*(-127140.0*exp(0.2444*V)-0.00003474*exp(-0.04391*V))+0.1212*exp(-0.01052*V)/(1.0+exp(-0.1378*(V+40.14))))

	E_K = R*T/F*log(K_o/K_i)
	i_t = g_to*r*(a*s+b*s_slow)*(V-E_K)

	r_inf = (1.0+exp((V+10.6)/(-11.42)))^(-1)
	tau_r = 100.0/(45.16*exp(0.03577*(V+50.0))+98.9*exp(-0.1*(V+38.0)))

	s_inf = (1.0+exp((V+45.3)/6.8841))^(-1)
	tau_s = 20.0*exp(-((V+70.0)/25.0)^2)+35.0

	s_slow_inf = (1.0+exp((V+45.3)/6.8841))^(-1)
	tau_s_slow = 1300.0*exp(-((V+70.0)/30.0)^2)+35.0

	i_ss = g_ss*r_ss*s_ss*(V-E_K)

	r_ss_inf = (1.0+exp((V+11.5)/(-11.82)))^(-1)
	tau_r_ss = 10000.0/(45.16*exp(0.03577*(V+50.0))+98.9*exp(-0.1*(V+38.0)))

	s_ss_inf = (1.0+exp((V+87.5)/10.3))^(-1)
	tau_s_ss = 2100.0

	i_K1 = 0.001*(0.048/(exp((V+37.0)/25.0)+exp((V+37.0)/(-25.0)))+0.01)/(1.0+exp((V-(E_K+76.77))/(-17.0)))+(g_K1*(V-(E_K+1.73)))/((1.0+exp(1.613*F*(V-(E_K+1.73))/(R*T)))*(1.0+exp((K_o-0.9988)/(-0.124))))

	f_K = 1.0-f_Na
	i_f_Na = g_f*y*f_Na*(V-E_Na)
	i_f_K = g_f*y*f_K*(V-E_K)
	i_f = i_f_Na+i_f_K

	y_inf = (1.0+exp((V+138.6)/10.48))^(-1)
	tau_y = 1000.0/(0.11885*exp((V+80.0)/28.37)+0.5623*exp((V+80.0)/(-14.19)))

	i_B_Na = g_B_Na*(V-E_Na)
	i_B_K = g_B_K*(V-E_K)

	nu_1 = gamma1*(1.0+KdNaio/Na_i)^2*(1.0+KdNais/Na_i*exp(-Delta*F*V/(R*T)))+(1.0+alpha/K_o)^2*(1.0+Na_o/KdNaes*exp(-(1.0-Delta)*F*V/(R*T)))
	i_NaK = i_NaK_max*(gamma1+1.0)/nu_1

	C_co = (Ca_i+J_R*Ca_SR/g_D)/(1.0+J_R/g_D)

	C_oc = (abs(FVRT_Ca) > 1e-9) ? (Ca_i+J_L/g_D*Ca_o*FVRT_Ca*exp(-FVRT_Ca)/(1.0-exp(-FVRT_Ca)))/(1.0+J_L/g_D*FVRT_Ca/(1.0-exp(-FVRT_Ca))) : ((Ca_i+J_L*Ca_o/g_D)/(1.0+J_L/g_D))

	expVL = exp((V-V_L)/(del_VL))

	alpha_p = expVL/(t_L*(1.0+expVL))
	alpha_m = phi_L/t_L

	beta_poc = C_oc^2/(t_R*(C_oc^2+K_RyR^2))
	beta_pcc = Ca_i^2/(t_R*(Ca_i^2+K_RyR^2))
	beta_m = phi_R/t_R

	epsilon_pco = C_co*(expVL+a1)/(tau_L*K_L*(1.0+expVL))
	epsilon_pcc = Ca_i*(expVL+a1)/(tau_L*K_L*(1.0+expVL))
	epsilon_m = b1*(expVL+a1)/(tau_L*(b1*expVL+a1))

	mu_poc = (C_oc^2+c*K_RyR^2)/(tau_R*(C_oc^2+K_RyR^2))
	mu_pcc = (Ca_i^2+c*K_RyR^2)/(tau_R*(Ca_i^2+K_RyR^2))
	mu_moc = theta_R*d*(C_oc^2+c*K_RyR^2)/(tau_R*(d*C_oc^2+c*K_RyR^2))
	mu_mcc = theta_R*d*(Ca_i^2+c*K_RyR^2)/(tau_R*(d*Ca_i^2+c*K_RyR^2))

	J_Rco = J_R*(Ca_SR-Ca_i)/(1.0+J_R/g_D)
	J_Roo = (abs(FVRT_Ca) > 1e-5) ? J_R*(Ca_SR-Ca_i+J_L/g_D*FVRT_Ca/(1.0-exp(-FVRT_Ca))*(Ca_SR-Ca_o*exp(-FVRT_Ca)))/(1.0+J_R/g_D+J_L/g_D*FVRT_Ca/(1.0-exp(-FVRT_Ca))) : J_R*(Ca_SR-Ca_i+J_L/g_D*1e-5/(1.0-exp(-1e-5))*(Ca_SR-Ca_o*exp(-1e-5)))/(1.0+J_R/g_D+J_L/g_D*1e-5/(1.0-exp(-1e-5)))
	J_Loc = (abs(FVRT_Ca) > 1e-5) ? J_L*FVRT_Ca/(1.0-exp(-FVRT_Ca))*(Ca_o*exp(-FVRT_Ca)-Ca_i)/(1.0+J_L/g_D*FVRT_Ca/(1.0-exp(-FVRT_Ca))) : J_L*1e-5/(1.0-exp(-1e-5))*(Ca_o*exp(-1e-5)-Ca_i)/(1.0+J_L/g_D*1e-5/(1.0-exp(-1e-5)))
	J_Loo = (abs(FVRT_Ca) > 1e-5) ? J_L*FVRT_Ca/(1.0-exp(-FVRT_Ca))*(Ca_o*exp(-FVRT_Ca)-Ca_i+J_R/g_D*(Ca_o*exp(-FVRT_Ca)-Ca_SR))/(1.0+J_R/g_D+J_L/g_D*FVRT_Ca/(1.0-exp(FVRT_Ca))) : J_L*1e-5/(1.0-exp(-1e-5))*(Ca_o*exp(-1e-5)-Ca_i+J_R/g_D*(Ca_o*exp(-1e-5)-Ca_SR))/(1.0+J_R/g_D+J_L/g_D*1e-5/(1.0-exp(-1e-5)))

	denom = (alpha_p+alpha_m)*((alpha_m+beta_m+beta_poc)*(beta_m+beta_pcc)+alpha_p*(beta_m+beta_poc))
	P1 = alpha_p*beta_m*(alpha_p+alpha_m+beta_m+beta_pcc)/denom
	P2 = alpha_m*(beta_pcc*(alpha_m+beta_m+beta_poc)+beta_poc*alpha_p)/denom
	P3 = alpha_p*(beta_poc*(alpha_p+beta_m+beta_pcc)+beta_pcc*alpha_m)/denom
	P4 = alpha_m*beta_m*(alpha_m+alpha_p+beta_m+beta_poc)/denom

	r_1 = P1*mu_poc+P4*mu_pcc
	r_2 = (alpha_p*mu_moc+alpha_m*mu_mcc)/(alpha_p+alpha_m)
	r_3 = beta_m*mu_pcc/(beta_m+beta_pcc)
	r_4 = mu_mcc
	r_5 = P2*epsilon_pco+P4*epsilon_pcc
	r_6 = epsilon_m
	r_7 = alpha_m*epsilon_pcc/(alpha_p+alpha_m)
	r_8 = epsilon_m

	z_4 = 1.0-z_1-z_2-z_3

	J_R1 = P3*J_Roo+J_Rco*P2
	J_R3 = J_Rco*beta_pcc/(beta_m+beta_pcc)
	i_RyR = N*(z_1*J_R1+z_3*J_R3)/Vmyo

	J_L1 = J_Loo*P3+J_Loc*P1
	J_L2 = J_Loc*alpha_p/(alpha_p+alpha_m)
	i_LCC1 = N*(z_1*J_L1+z_2*J_L2)/Vmyo
	i_LCC2 = -i_LCC1*2.0*F*Vmyo_uL

	i_NCX1 = g_NCX*(Na_i^3*Ca_o*exp(eta*FVRT)-Na_o^3*Ca_i*exp((eta-1.0)*FVRT))/((Na_o^3+K_mNa^3)*(Ca_o+K_mCa)*(1.0+k_sat*exp((eta-1.0)*FVRT)))
	i_NCX2 = i_NCX1*F*Vmyo_uL

	i_SERCA = g_SERCA*Ca_i^2/(K_SERCA^2+Ca_i^2)

	i_pCa1 = g_pCa*Ca_i/(K_mpCa+Ca_i)
	i_pCa2 = i_pCa1*2*F*Vmyo_uL

	E_Ca = R*T*log(Ca_o/Ca_i)/(2.0*F)
	i_CaB1 = g_CaB*(E_Ca-V)
	i_CaB2 = -i_CaB1*2.0*F*Vmyo_uL

	i_SR = g_SRl*(Ca_SR-Ca_i)

	i_TRPN = k_m_TRPN*(B_TRPN-TRPN)-k_p_TRPN*TRPN*Ca_i

	beta_CMDN = (1.0+k_CMDN*B_CMDN/(k_CMDN+Ca_i)^2+EGTA_tot*KmEGTA/(KmEGTA+Ca_i)^2)^(-1)

	du[1] = dz_1 = -(r_1 + r_5)*z_1 + r_2*z_2 + r_6*z_3
	du[2] = dz_2 = r_1*z_1 - (r_2 + r_7)*z_2 + r_8*z_4
	du[3] = dz_3 = r_5*z_1 - (r_6 + r_3)*z_3 + r_4*z_4
	du[4] = dr = (r_inf - r)/tau_r
	du[5] = ds = (s_inf - s)/tau_s
	du[6] = ds_slow = (s_slow_inf - s_slow)/tau_s_slow
	du[7] = dy = (y_inf - y)/tau_y
	du[8] = dCa_SR = Vmyo_uL*(-i_RyR + i_SERCA - i_SR)/VSR_uL
	du[9] = dCa_i = beta_CMDN*(i_RyR - i_SERCA + i_SR + i_TRPN - (-2.0*i_NCX2 + i_LCC2 + i_pCa2 + i_CaB2)/(2.0*Vmyo_uL*F))
	du[10] = dK_i = -(i_Stim + i_ss + i_B_K + i_t + i_K1 + i_f_K -2.0*i_NaK)/(Vmyo_uL*F)
	du[11] = dNa_i = -(i_Na + i_B_Na + 3.0*i_NCX2 + 3.0*i_NaK + i_f_Na)/(Vmyo_uL*F)
	du[12] = dTRPN = i_TRPN
	du[13] = dV = -(i_Na + i_t + i_ss + i_f + i_K1 + i_B_Na + i_B_K + i_NaK + i_CaB2 + i_NCX2 + i_pCa2 + i_LCC2 + i_Stim)/Cm
	du[14] = dh = (h_inf - h)/tau_h
	du[15] = dj = (j_inf - j)/tau_j
	du[16] = dm = (m_inf - m)/tau_m
	du[17] = dr_ss = (r_ss_inf - r_ss)/tau_r_ss
	du[18] = ds_ss = (s_ss_inf - s_ss)/tau_s_ss
end

# Solve the ODEs system
function solve_ep(rat, freq, n_beats, n_curves, soltype, args...)
	p = init_constants(rat, freq)
	stim_period = p[5]

	# Available additional argument is the vector of pertubations
	if length(args) > 0
		S = args[1]
		l = [11, 13, 17, 18, 19, 23, 34, 54, 55, 57, 60] # indices of proteins to perturbate
		p[l] = S.*p[l] # perform the actual perturbation by means of a scaling vector S
	end
	
	tspan = (0.0, n_beats*stim_period-1.0)

	dt = soltype == "fragmentary" ? stim_period : 1.0
	t = collect((n_beats-n_curves)*stim_period:dt:n_beats*stim_period-1.0) # points where to save the sol.

	u0 = init_states(rat, freq)
	prob = ODEProblem(fun, u0, tspan, p)
	
	alg = CVODE_BDF()
	sol = solve(prob, alg, dtmax=1.0, saveat=t)

	return t.-(n_beats-n_curves)*stim_period, sol
end

end # end module