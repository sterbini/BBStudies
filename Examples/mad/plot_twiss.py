
import matplotlib.pyplot as plt
def pltwis(tw):
	# TWISS PLOT BLOCK
	fig1 = plt.figure(1, figsize=(6.4, 4.8*1.5))
	spbet = plt.subplot(3,1,1)
	spco = plt.subplot(3,1,2, sharex=spbet)
	spdisp = plt.subplot(3,1,3, sharex=spbet)

	spbet.plot(tw['s'], tw['betx'])
	spbet.plot(tw['s'], tw['bety'])

	spco.plot(tw['s'], tw['x'])
	spco.plot(tw['s'], tw['y'])

	spdisp.plot(tw['s'], tw['dx'])
	spdisp.plot(tw['s'], tw['dy'])


	spbet.set_ylabel(r'$\beta_{x,y}$ [m]')
	spco.set_ylabel(r'(orbit)$_{x,y}$ [m]')
	spdisp.set_ylabel(r'$D_{x,y}$ [m]')
	spdisp.set_xlabel('s [m]')
	#fig1.subplots_adjust(left=.15, right=.92, hspace=.27)
