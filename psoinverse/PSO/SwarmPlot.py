import matplotlib.pyplot as pl
import matplotlib.animation as man
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

class SwarmPlotter(object):
   # TODO: add path, rerun LAM with new structure of evaluate, and figure out how to keep a record of the PBest run.
    def __init__(self, runname, numplotcoords, maxnsteps, labels=["x","y","z"], axislimits=[[0,1],[0,1],[0,1]], OptimalCoord=None, makemovie=False, OutputDir=None, **kwargs):
        self.numplotcoords = numplotcoords
        if numplotcoords > 3:
            print "Bad numplotcoords = {}. Only first three coordinates will be plotted".format(numplotcoords)
            self.numplotcoords = 3
        elif numplotcoords < 1:
            print "Bad numplotcoords = {}. First coordinate will be plotted.".format(numplotcoords)
            self.numplotcoords = 1
        self.labels = labels
        self.runname = runname
        self.OutputDir = OutputDir
        self.OptimalCoord = OptimalCoord

        # Initialize the plot
        self.fig = pl.figure(figsize=(8,4))
        if self.numplotcoords == 3:
            self.ax_left = self.fig.add_subplot(121,projection='3d')
        else:
            self.ax_left = self.fig.add_subplot(121)
        self.ax_right = self.fig.add_subplot(122)

        # Set axis labels and scale
        self.axislimits = axislimits

        self.ax_right.set_xlabel("Iteration")
        self.ax_right.set_xlim(0,maxnsteps)
        self.ax_right.set_ylabel("Fitness")


        self.makemovie = makemovie
        if makemovie:
            FFMPEGwriter = man.writers['ffmpeg']
            metadata = dict(title=self.runname)
            fps = kwargs.get("fps",8)
            dpi = kwargs.get("dpi",100)
            self.writer = FFMPEGwriter(fps=fps, metadata=metadata)
            if self.OutputDir == None:
                self.writer.setup(self.fig, "{}.mp4".format(self.runname), dpi=dpi)
            else:
                self.writer.setup(self.fig, "{}/{}.mp4".format(self.OutputDir,self.runname), dpi=dpi)

    def updatePlot(self, agent_list, gbest, currentstep, overlaycallbackfn=None, fitnessplotmin=-1e2, **kwargs):
        # The left panel depends on dimensionality
        if self.numplotcoords == 1:
            self.__updatePlot1D(agent_list, gbest, currentstep, overlaycallbackfn, **kwargs)
        elif self.numplotcoords == 2:
            self.__updatePlot2D(agent_list, gbest, currentstep, overlaycallbackfn, **kwargs)
        elif self.numplotcoords == 3:
            self.__updatePlot3D(agent_list, gbest, currentstep, overlaycallbackfn, **kwargs)
        else:
            print "Unknown plot type"

        # The right panel is always the same - the fitness over time
        # BUT here we mask out extreme fitness values that are sometimes used for
        #  invalid / pathological evaluations in computeFitness
        for a in agent_list:
            if a.Location.Fitness > fitnessplotmin:
                self.ax_right.scatter(currentstep, a.Location.Fitness, c='grey', alpha=0.4)
        if gbest.Fitness > fitnessplotmin:
            self.ax_right.scatter(currentstep, gbest.Fitness, c='b')

        # Compact the layout and update the plot and, if enabled, the movie
        self.fig.tight_layout()
        if self.OutputDir == None:
            self.fig.savefig("{}.pdf".format(self.runname))
        else:
            self.fig.savefig("{}/{}.pdf".format(self.OutputDir,self.runname))
        if self.makemovie:
            self.writer.grab_frame()
        pl.show(block=False)


    def __updatePlot1D(self, agent_list, gbest, currentstep, overlaycallbackfn=None, **kwargs):
        self.ax_left.clear()
        self.ax_left.set_xlabel("Agent ID")
        self.ax_left.set_xlim(0.5,len(agent_list)+0.5)
        self.ax_left.set_ylabel(self.labels[0])
        self.ax_left.set_ylim(self.axislimits[0],self.axislimits[1])

        # The caller can provide a callback function to add something to the axes
        if overlaycallbackfn != None:
            overlaycallbackfn(self.ax_left, **kwargs)

        # Reset the series cycler
        self.ax_left.set_prop_cycle(None)
        j=1
        for a in agent_list:
            x = a.get_coords()[0]
            vx = a.Velocity[0]*a.Location.Scale[0]
            self.ax_left.scatter(j,x)
            self.ax_left.arrow(j,x,0,vx,head_width = 0.05)
            j+=1
        # Reset the series cycler
        self.ax_left.set_prop_cycle(None)
        j=1
        for a in agent_list:
            x = a.PBest.get_scaled_coords()[0]
            self.ax_left.scatter(j,x,alpha=0.2) # Faded point with no v arrow for pbest
            j+=1

        # Draw a line for the gbest solution
        gbestval = gbest.get_scaled_coords()[0]
        self.ax_left.plot([0,len(agent_list)+1],[gbestval,gbestval],c='#222222',ls='-',linewidth=1.,label='gbest (to date)')
        if self.OptimalCoord != None:
            self.ax_left.plot([0,len(agent_list)+1],[self.OptimalCoord,self.OptimalCoord],c='#D4AF37',ls='-',linewidth=2.,label='SCFT varcell')
        self.ax_left.legend(loc=1)


    def __updatePlot2D(self, agent_list, gbest, currentstep, overlaycallbackfn=None, **kwargs):
        self.ax_left.clear()
        self.ax_left.set_xlabel(self.labels[0])
        self.ax_left.set_xlim(self.axislimits[0][0],self.axislimits[0][1])
        self.ax_left.set_ylabel(self.labels[1])
        self.ax_left.set_ylim(self.axislimits[1][0],self.axislimits[1][1])

        # The caller can provide a callback function to add something to the axes
        if overlaycallbackfn != None:
            overlaycallbackfn(self.ax_left, **kwargs)

        # Reset the series cycler
        self.ax_left.set_prop_cycle(None)
        for a in agent_list:
            x = a.get_coords()[0]
            y = a.get_coords()[1]
            vx = a.Velocity[0]*a.Location.Scale[0]
            vy = a.Velocity[1]*a.Location.Scale[1]
            self.ax_left.scatter(x,y)
            # ax.arrow is affected by axis limits and aspect ratio. Therefore we use annotate (as recommended by documentation)
            # self.ax_left.arrow(x,y,vx,vy,width=0.1)
            AN = self.ax_left.annotate("", xy=(x+vx,y+vy), xycoords='data', xytext=(x,y), textcoords='data', arrowprops=dict(arrowstyle="-|>", fc="k"), annotation_clip=False)
            AN.arrow_patch.set_clip_box(self.ax_left.bbox)
        # Reset the series cycler
        self.ax_left.set_prop_cycle(None)
        for a in agent_list:
            x = a.PBest.get_scaled_coords()[0]
            y = a.PBest.get_scaled_coords()[1]
            self.ax_left.scatter(x,y,alpha=0.2) # Faded point with no v arrow for pbest

        # Put a ring around the fittest agent
        self.ax_left.scatter(gbest.get_scaled_coords()[0],gbest.get_scaled_coords()[1],s=300,edgecolors='r',c='',linewidths=10.)

        # Render the target solution, if available
        if self.OptimalCoord != None:
            self.ax_left.scatter(self.OptimalCoord[0],self.OptimalCoord[1],s=100,marker='*',c='#D4AF37')

    def __updatePlot3D(self, agent_list, gbest, currentstep, overlaycallbackfn=None, **kwargs):
        self.ax_left.clear()
        self.ax_left.set_xlabel(self.labels[0])
        self.ax_left.set_xlim(self.axislimits[0][0],self.axislimits[0][1])
        self.ax_left.set_ylabel(self.labels[1])
        self.ax_left.set_ylim(self.axislimits[1][0],self.axislimits[1][1])
        self.ax_left.set_zlabel(self.labels[2])
        self.ax_left.set_zlim(self.axislimits[2][0],self.axislimits[2][1])

        # The caller can provide a callback function to add something to the axes
        if overlaycallbackfn != None:
            overlaycallbackfn(self.ax_left, **kwargs)

        # Optionally could color by fitness:
        # vmin = kwargs.get("vmin", np.min(c))
        # vmax = kwargs.get("vmax", np.max(c))
        # c_scaled = (c-vmin)/(vmax-vmin)
        # Loop points xv, yv, zv, cv in x,y,z,c_scaled
        #   alpha = max(.1, min(1, 1-cv))
        #   axis.scatter(xv, yv, zv, c=cv, alpha=alpha, norm=mpl.colors.Normalize(vmin=0, vmax=1), depthshade=False)

        # Reset the series cycler
        self.ax_left.set_prop_cycle(None)
        j=1
        for a in agent_list:
            x = a.get_coords()[0]
            y = a.get_coords()[1]
            z = a.get_coords()[2]
            vx = a.Velocity[0] * a.Location.Scale[0]
            vy = a.Velocity[1] * a.Location.Scale[1]
            vz = a.Velocity[2] * a.Location.Scale[2]
            self.ax_left.scatter(x,y,z)
            self.ax_left.quiver(x,y,z,vx,vy,vz)
            j+=1
        # Reset the series cycler
        self.ax_left.set_prop_cycle(None)
        j=1
        for a in agent_list:
            x = a.PBest.get_scaled_coords()[0]
            y = a.PBest.get_scaled_coords()[1]
            z = a.PBest.get_scaled_coords()[2]
            self.ax_left.scatter(x,y,z,alpha=0.2) # Faded point with no v arrow for pbest
            j+=1

        # Draw a line for the gbest solution
        gbestval = gbest.get_scaled_coords()
        self.ax_left.scatter(gbestval[0],gbestval[1],gbestval[2],s=300,edgecolors='r',c='',linewidths=10.)
        if self.OptimalCoord != None:
            self.ax_left.scatter(self.OptimalCoord[0],self.OptimalCoord[1],self.OptimalCoord[2],s=100,marker='*',c='#D4AF37')


    def finishMovie(self):
        if self.makemovie:
            self.writer.finish()



