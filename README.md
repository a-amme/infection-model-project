# Replication and Optimization of a Multicellular Spatial Model of IFN Response to Viral Infection

This respository holds materials for our replication and optimization of the simulation reported in the paper
**Multicellular Spatial Model of RNA Virus Replication and Interferon Responses Reveals Factors Controlling Plaque Growth Dynamics**: 
https://www.biorxiv.org/content/10.1101/2021.03.16.435618v1

The original model code can be found here on GitHub in the repository **Multicellular Spatial Model of RNA Virus Replication**:
https://github.com/ImmuSystems-Lab/Multicellular_Spatial_Model_of_RNA_Virus_Replication

This model is designed to be accessible regardless of experience with computer programming or CompuCell3D. To facilitate exploration for the curious, we suggest the following exercises: 

<ol>
  <li>
    In this model, diffusion of virions is the only way the virus can move around the 
    host cell culture. The tendency of the virions to diffuse is modeled mathematically 
    by the viral diffusion constant, which you can change while the simulation runs 
    using the sliders. How would the infection progress if the virus diffused more than 
    the default? Or less? Use the sliders to change the value of the viral diffuison 
    constant – as one would expect, a higher viral diffusion constant means more diffusion, 
    while a lower constant means less diffusion – and compare. 
  </li>
  <li>
    Virions don't just move through the extracellular matrix by diffusion in this model; 
    they are also uptaken by host cells, which we capture mathematically as virion decay.
    What if host cells were to take up greater concentrations of virions at once, or 
    something were to destroy the virions in the extracellular matrix? You can experiment 
    with this situation by using the sliders to intensify the viral decay. Try increasing 
    viral decay by a factor of 10. How does the simulation's behavior change? Why? 
  </li>
  <li>
    Once simulated cells enter the viral-releasing stage, the likelihood that they will 
    transition to the death state increases as their health decreases. What would happen 
    if the virus were more virulent than the default? Increase the rate at which infected cells 
    lose health and compare what you see to the default simulation's behavior. Does the 
    infection spread more or less? Why? 
  </li>
  <li>
    The cells in this model alert one another to the presence of the virus by releasing 
    interferons. The interferon-releasing process within the cell is accelerated when the 
    cell detects extracellular interferon, and this detection is accomplished by the 
    molecule STAT, which "activates" and becomes STATP upon interaction with interferon. 
    Increase STAT activation by a factor of 100 and compare the infection outcomes. What 
    is the difference? Why? 
  </li>
  <li>
    The cells in this model alert one another to the presence of the virus by releasing 
    interferons. What if the interferon signal were more potent, spreading throughout 
    more of the simulated cell culture? Increase the diffuison constant of extracellular 
    interferon using the steering panels and compare the simulation's behavior with that 
    of the default. 
  </li>
  <li>
    Researchers are working on treatments in which infected tissues are bathed in 
    cytokines such as interferons prior to or at the onset of infection. It's possible 
    to simulate such a treatment with this model: open the model in CompuCell3D's 
    Twedit code editor, find the varible <code>IFNWash</code>, and set it to <code>True</code>. 
    How does the infection progress under these conditions? Based on your observations 
    alone, is this a good treatment? Would you expect this treatment to be as 
    effective or as ineffective as this if it were done <i>in vivo</i>? Why or why not? 
  </li>
</ol>
