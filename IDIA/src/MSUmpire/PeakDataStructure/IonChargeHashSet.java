/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package MSUmpire.PeakDataStructure;

/**
 *
 * @author
 */
public final class IonChargeHashSet {
    private short state=0;
    public void add(final int charge){
        assert charge>0 && charge<16 : "charge must be in (0,16):\t" + charge;
        this.state |= 1 << charge;
    }
    public boolean contains(final int charge){
        assert charge>0 && charge<16 : "charge must be in (0,16):\t" + charge;
//        return (this.state & 1 << charge) != 0;
        return (this.state >>> charge & 1) != 0;
    }

    public static short add(final short state,final int charge){
        assert charge>0 && charge<16 : "charge must be in (0,16):\t" + charge;
        return (short) (state | 1 << charge);
    }
    public static boolean contains(final short state, final int charge){
        assert charge>0 && charge<16 : "charge must be in (0,16):\t" + charge;
        return (state >>> charge & 1) != 0;
    }

}
