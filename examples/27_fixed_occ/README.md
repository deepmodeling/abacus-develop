For a spin-polarized calculation (where `nspin=2` in the INPUT), let's say you want to excite a spin-up (spin component 1) or spin-down (spin component 2) electron from the ground state to the excited state.

To do this, you will need to do the following:

- Find out how many bands do you have. You can obtain that by search for `NBANDS` in the running_*.log file
- Find out how many k points do you have.  Let's call that number `N_K`.
- Find out how many electrons are occupying spin components 1 and 2. Let's say there are `O_UP` electrons that occupy the spin up, `O_DN` electrons that occupy the spin down.
- Next, calculate how many states are empty in both spins: empty states in the spin up = `U_UP` where `U_UP = NBANDS - O_UP`. Same applies for `U_DN`.

Finally, add the `ocp_set` as follows:

- `ocp_set  [O_UP-1]*1.0 1*0.0 1*1.0 [U_UP-1]*0.0 <repeated N_K times> [O_DN]*1.0 [U_DN]*0.0 <repeated N_K times>`

Note the `<repeated N_K times>` above. It means: repeat `[O_UP-1]*1.0 1*0.0 1*1.0 [U_UP-1]*0.0` for `N_K` times.

Let's have an example. Let's say we have 8 kpoints in both spins, 216 bands, where the spin up electrons occupy 144 bands and the spin down occupy 142 bands. Then, here are the two tags:

- `ocp_set 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 143*1.0 1*0.0 1*1.0 71*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0 142*1.0 74*0.0`

I hope that helps!
