# ViralDrive
Viral gene drive simulation

This is a R program that perform simple simulation of a viral gene drive. See paper (once it is submitted).
This is really a work in progress, insuffiently commented i'm aware. So please ask questions and i'll do my best to help.

In brief, the simulation of works as follow.
4 type of viruses exist:
- WT: wildtype virus that can be converted to new gene drive virus
- GD: Gene drive virus. If coinfecting a cell with a WT virus, it convert the WT virus into GDo or Res virus
- GDo: Same as GD essentially, except that it is marked as originating from a WT virus. It can convert new WT virus
- R: Virus resistant to the drive, don't care about GD and GDo

GD and GDo replicate with a fitness cost f, and convert WT viruses with an efficiency 1-alpha

At each viral generation, NxMOI virus infect N cells
- If a cell is infected by WT and/or R virus, it produce in average 100 (SD 30, gaussian) virus in the same proportion
- If a cell is infected by GD and/or GDo virus, it produce in average fx100 virus in the same proportion. 
- If a cell is infected by WT AND GD or GDo, 100 virus are created, but WT viruses are converted into GDo (proportion 1-alpha) and R virus (proportion alpha)

After doing that for every cell at a given generation, we calculate the new proportion of the different viruses, and use it to iter to the next generation
