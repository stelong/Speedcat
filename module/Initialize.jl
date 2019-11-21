# Define module
module Initialize

# Make them visible outside
export init_constants
export init_states

function init_constants(rat, freq)
	if rat != "sham" && rat != "ab"
		println("\nInvalid rat phenotype! Try with either \"sham\" or \"ab\".\n")
		return 0
	end
	if freq != 1 && freq != 6
		println("\nInvalid frequency! Try with either \"1\" or \"6\".\n")
		return 0
	end

	R               = 8314.0
	T               = 310.0
	F               = 96487.0
	Cm              = 0.0001
	stim_duration   = 3.0
	stim_amplitude  = -0.0012
	Vmyo            = 25850.0
	Vmyo_uL         = 2.59e-5
	VSR_uL          = 2.10e-6
	Na_o            = 140.0
	a               = 0.883
	b               = 0.117
	K_o             = 5.4
	g_f             = 1.45e-6
	f_Na            = 0.2
	g_B_Na          = 8.02e-8
	g_B_K           = 1.38e-7
	i_NaK_max       = 0.00138
	gamma1          = 3.6
	KdNaio          = 19.0
	KdNais          = 22.0
	KdNaes          = 880.0
	Delta           = 0.3
	alpha           = 1.8
	Ca_o            = 1.8
	g_D             = 0.1
	J_R             = 0.02
	N               = 50000.0
	del_VL          = 7.0
	phi_L           = 11.5
	t_L             = 1.0
	tau_L           = 1550.0
	t_R             = 1.17
	tau_R           = 2.4
	phi_R           = 0.05
	theta_R         = 0.012
	K_RyR           = 0.065
	a1              = 0.0625
	b1              = 14.0
	c               = 0.01
	d               = 100.0
	K_mNa           = 87.5
	K_mCa           = 1.38
	eta             = 0.35
	k_sat           = 0.1
	g_pCa           = 5.00e-6
	K_mpCa          = 0.00035
	g_SRl           = 1.00e-6
	k_m_TRPN        = 0.04
	k_p_TRPN        = 40.0
	B_TRPN          = 0.07
	k_CMDN          = 0.002382
	B_CMDN          = 0.05
	EGTA_tot        = 0.0
	KmEGTA          = 0.00015
	tau_s_ss        = 2100.0
	f_K             = 1.00000 - f_Na
	alpha_m         = phi_L/t_L
	beta_m          = phi_R/t_R

	if rat == "sham"
		g_Na = 0.0007
		g_to = 2.0e-5
		g_ss = 1.3e-5
		g_K1 = 4.0e-5
		J_L = 0.0008
		V_L = -9.0
		K_L = 0.00038
		g_NCX = 0.0234
		g_CaB = 2.0e-8

		if freq == 1
			stim_period = 1000.0
			g_SERCA = 0.000235
			K_SERCA = 0.0004968
		else
			stim_period = 170.0
			g_SERCA = 0.00051
			K_SERCA = 0.00069
		end
	else
		g_Na = 0.0002
		g_to = 1.4e-5
		g_ss = 1.0e-6
		g_K1 = 1.5e-5
		J_L = 0.0012 
		V_L = -13.0
		K_L = 0.00016
		g_NCX = 0.0456
		g_SERCA = 0.00049
		g_CaB = 6.0e-9

		if freq == 1
			stim_period = 1000.0
			K_SERCA = 0.00044
		else
			stim_period = 170.0
			K_SERCA = 0.00025
		end
	end

	p = [R, T, F, Cm, stim_period, stim_duration, stim_amplitude, Vmyo, Vmyo_uL, VSR_uL, g_Na, Na_o, g_to, a, b, K_o, g_ss, g_K1, g_f, f_Na, g_B_Na, g_B_K, i_NaK_max, gamma1, KdNaio, KdNais, KdNaes, Delta, alpha, Ca_o, g_D, J_R, J_L, N, V_L, del_VL, phi_L, t_L, tau_L, t_R, tau_R, phi_R, theta_R, K_RyR, K_L, a1, b1, c, d, K_mNa, K_mCa, eta, k_sat, g_NCX, g_SERCA, K_SERCA, g_pCa, K_mpCa, g_CaB, g_SRl, k_m_TRPN, k_p_TRPN, B_TRPN, k_CMDN, B_CMDN, EGTA_tot, KmEGTA, tau_s_ss, f_K, alpha_m, beta_m]
	return p
end


function init_states(rat, freq)
	if rat != "sham" && rat != "ab"
		println("\nInvalid rat phenotype! Try with either \"sham\" or \"ab\".\n")
		return 0
	end
	if freq != 1 && freq != 6
		println("\nInvalid frequency! Try with either \"1\" or \"6\".\n")
		return 0
	end

	if rat == "sham"
		if freq == 1
			u0 = [0.990016532916529, 0.00845823628523856, 0.00151233172289407, 0.00146331830465093, 0.996934138278418, 0.78841193673441, 0.00566123148325894, 0.96268028201207, 0.000103020385969363, 142.919492013701, 8.46899983583716, 0.0633670056927004, -85.1221681609219, 0.815520320018128, 0.815471795073686, 0.00208137744708665, 0.0019683140031203, 0.416987850222633]
		else
			u0 = [0.915759670501361, 0.0116308821918943, 0.0716996181298885, 0.00140717973440432, 0.954569649396715, 0.300568843000913, 0.00395453931879583, 1.0663135643117, 0.000463188750506148, 139.994789812443, 12.960204722764, 0.0462772817360074, -85.5696059816548, 0.824266840236936, 0.722995966352439, 0.00194320513549252, 0.00189620611861644, 0.310616913661423]
		end
	else
		if freq == 1
			u0 = [0.989618266628688, 0.00828851292530188, 0.00207588449544264, 0.00150276538930736, 0.996791125715503, 0.625815901322428, 0.00484325258805042, 1.22566116672694, 4.53106928940201e-5, 139.964416466575, 7.86198765026574, 0.0669207270187171, -84.8179335690707, 0.807799685437556, 0.807494936350446, 0.00218089269066671, 0.0020195578095622, 0.368560872566697]
		else
			u0 = [0.650490049579104, 0.00570189155987283, 0.340821224038804, 0.00132309740749641, 0.913059928712858, 0.215036691884835, 0.0025539544569868, 1.78499523287115, 0.000145080891755077, 140.1638657529, 10.83563050735, 0.0590608556435992, -86.2742001770196, 0.831855538514125, 0.585528587217805, 0.00174392847776228, 0.00178817500396492, 0.209855226174927]
		end
	end
	return u0
end

end # end module