#include "Orthogonal.h"
#include "tools.h"

Orthogonal::Orthogonal()
{
    test=0;
}
Orthogonal::~Orthogonal() {}

void Orthogonal::start(SpillageStep &step1, SpillageStep &step2)
{
    if (test==1)TITLE("Orthogonal","start");
    timer::tick("Orthogonal","start");

    Inverse_Matrix inverse;
    for (int istr=0; istr<STRNUM; istr++)
    {
        assert( step1.data[istr].nks == step2.data[istr].nks );
        assert( step1.data[istr].nbands == step2.data[istr].nbands );
        assert( step1.data[istr].ne == step2.data[istr].ne );
        assert( step1.data[istr].nwfc_all == step2.data[istr].nwfc_all );

        const int nks = step1.data[istr].nks;
        const int nbands = step1.data[istr].nbands;
        const int ne = step1.data[istr].ne;
        const int nwfc_all = step1.data[istr].nwfc_all;
        const int nwfc2 = step1.data[istr].nwfc2;

        if (test==1)
        {
            cout << "\n" << setw(10) << "nks"
                 << setw(10) << "nbands"
                 << setw(10) << "ne"
                 << setw(10) << "nwfc_all"
                 << setw(10) << "nwfc2";

            cout << "\n" << setw(10) << nks
                 << setw(10) << nbands
                 << setw(10) << ne
                 << setw(10) << nwfc_all
                 << setw(10) << nwfc2;
        }

        inverse.init( step1.data[istr].nwfc2 );
        for (int ik=0; ik<nks; ik++)
        {
//			if(ik==0)
//			{
//				PRINTCM("Sinv",step1.data[istr].Sinv[ik]);
//			}
            inverse.using_zheev( step1.data[istr].Soverlap[ik], step1.data[istr].Sinv[ik]);
//			if(ik==0)
//			{
//				PRINTCM("Sinv",step1.data[istr].Sinv[ik]);
//			}

                complex<double> *sum = new complex<double>[nbands];
                for (int iw=0; iw<nwfc_all; iw++)
                {
                    for (int ie=0; ie<ne; ie++)
                    {
                        ZEROS(sum, nbands);
                        for (int mu=0; mu<nwfc2; mu++)
                        {
                            const int Tmu = step1.wayd[mu].type;
                            const int Lmu = step1.wayd[mu].L;
                            const int Nmu = step1.wayd[mu].N;
                            const int iw00 = step1.wayd[mu].iw00;

                            complex<double> first_part = complex<double>(0,0);
							if(USEPW)
							{
								const int Mmu = step1.wayd[mu].m;
								const int Imu = step1.wayd[mu].i;
								PW.calculate_Jlq(ik,iw,ie);
								first_part = PW.calculate_Jlq_Phi(ik,mu);
								
						//		if( norm(first_part) > 0.1 )
						//		{
						//			cout << "\n iw=" << iw << " ie=" << ie << " mu=" << mu << " first_part=" << first_part;
						//			BLOCK_HERE("haha");
						//		}
							}
							else
							{
                            	for (int iemu=0; iemu<ne; iemu++)
                            	{
                                	first_part += input.Coef.C4( Tmu, Lmu, Nmu, iemu )
                                              * input.QS_data[istr].Sq1q2[ik]( iw, iw00, ie, iemu );
								}
                            }

                            for (int nu=0; nu<nwfc2; nu++)
                            {
                                for (int ib=0; ib<nbands; ib++)
                                {
                                    sum[ib] += first_part * step1.data[istr].Sinv[ik](mu,nu) *
                                               step1.data[istr].Qoverlap[ik](nu,ib);
                                }
                            }
                        }

                        for (int ib=0; ib<nbands; ib++)
                        {
                            step2.data[istr].Qin(ik,ib,iw,ie) =
                                step1.data[istr].Qin(ik,ib,iw,ie) - sum[ib];
                        }

                    }
                }
                delete[] sum;

			/*


            // old version : band loop should not be outside.
            for(int ib=0; ib<nbands; ib++)
            {
            	for(int iw=0; iw<nwfc_all; iw++)
            	{
            		for(int ie=0; ie<ne; ie++)
            		{
            			complex<double> sum = complex<double>(0,0);
            			for(int mu=0; mu<nwfc2; mu++)
            			{
            				const int Tmu = step1.wayd[mu].type;
            				const int Lmu = step1.wayd[mu].L;
            				const int Nmu = step1.wayd[mu].N;
            				const int iw00 = step1.wayd[mu].iw00;

            				for(int nu=0; nu<nwfc2; nu++)
            				{
            					complex<double> first_part = complex<double>(0,0);

            					for(int iemu=0; iemu<ne; iemu++)
            					{
            						first_part += input.Coef.C4( Tmu, Lmu, Nmu, iemu )
            							* input.QS_data[istr].Sq1q2[ik]( iw, iw00, ie, iemu );
            					}
            					sum += first_part * step1.data[istr].Sinv[ik](mu,nu) *
            						step1.data[istr].Qoverlap[ik](nu,ib);
            				}
            			}
            			step2.data[istr].Qin(ik,ib,iw,ie) =
            				step1.data[istr].Qin(ik,ib,iw,ie) - sum;
            		}// end ie
            	}// end iw
            }// end ib

			*/
        }
    }

    timer::tick("Orthogonal","start");
    return;
}
