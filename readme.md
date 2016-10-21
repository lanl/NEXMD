  `.. ..      `....    `..       `..`...     `..      `.       `........  `.. ..  `..       `..`.....    
`..    `..  `..    `.. `. `..   `...`. `..   `..     `. ..     `..      `..    `..`. `..   `...`..   `.. 
 `..      `..       `..`.. `.. ` `..`.. `..  `..    `.  `..    `..       `..      `.. `.. ` `..`..    `..
   `..    `..       `..`..  `..  `..`..  `.. `..   `..   `..   `......     `..    `..  `..  `..`..    `..
      `.. `..       `..`..   `.  `..`..   `. `..  `...... `..  `..            `.. `..   `.  `..`..    `..
`..    `..  `.. `. `.. `..       `..`..    `. .. `..       `.. `..      `..    `..`..       `..`..   `.. 
  `.. ..      `.. ..   `..       `..`..      `..`..         `..`........  `.. ..  `..       `..`.....    
                   `.                                                                                    
                   
This is the master version of SQM-NAESMD beta 1.0. See the manual for description, compilation, and operation of the program. You will need to have a compiled version of BLAS and LAPACK libraries to compile this code, which are likely already on your system and linked to by your compiler. Otherwise, get it at http://www.netlib.org/lapack/.

Bug list: Nonesofar

Several development branches of this program exist. These include:

amber_integration - working towards an integration with Amber, currently bootstrapping to Amber for NAESMD-QM/MM MD
correlation - attempting to develop a self-consistent TD-SCF method
eigensolvers - reorganizing the code to accept alternatives to the build in Lancsoz algorithm and Liouville operator sparsity
esxlbomd - development of excited state extended Lagrangian molecular dynamics
planning - a junk branch for testing random ideas
release0 - historical branch 1
release1 - historical branch 2
ssdynamics - nonequilibrium solvation models for molecular dynamics as published
xlbomd - development of ground state extended Lagrangian molecular dynamics

Josiah Bjorgaard, 2016, LANL
