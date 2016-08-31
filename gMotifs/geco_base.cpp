/***************************************************************************
 *   Copyright (C) 2010 by Uberto Pozzoli and Matteo Cereda                *
 *   uberto.pozzoli@bp.lnf.it                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "geco_base.h"
//#include <string.h>
#include <cstring>

namespace geco {
//--------------------------------------------------------------------------
// helper bits function definition
//--------------------------------------------------------------------------
const _unsigned_char mLogNBits[64]={
    0,0,0,0,0,0,0,3,
    0,0,0,0,0,0,0,4,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,5,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,6
};
const _unsigned_char mNBytes=sizeof(gBytes);
const _unsigned_char mNBits=mNBytes<<3;
const _unsigned_char mP1=mNBits-1;
const _unsigned_char mP2= mLogNBits[mP1];
const gBytes mP3= (gBytes) 1 << mP1;
#ifdef __GECO_GBYTES_UCHAR__
const gBytes MG[8]={128,64,32,16,8,4,2,1};
const gBytes MR[8]={127,191,223,239,247,251,253,254};
#endif
#ifdef __GECO_GBYTES_USHORT__
const gBytes MG[16]={32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1};
const gBytes MR[16]={32767,49151,57343,61439,63487,64511,65023,65279,65407,65471,65503,65519,65527,65531,65533,65534};
#endif
#ifdef __GECO_GBYTES_UINT__
const gBytes MG[32]={2147483648,1073741824,536870912,268435456,134217728,67108864,33554432,16777216,8388608,4194304,2097152,1048576,524288,262144,131072,65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1};
const gBytes MR[32]={2147483647,3221225471,3758096383,4026531839,4160749567,4227858431,4261412863,4278190079,4286578687,4290772991,4292870143,4293918719,4294443007,4294705151,4294836223,4294901759,4294934527,4294950911,4294959103,4294963199,4294965247,4294966271,4294966783,4294967039,4294967167,4294967231,4294967263,4294967279,4294967287,4294967291,4294967293,4294967294};
#endif
#ifdef __GECO_GBYTES_ULONG__
const gBytes MG[64]={0x8000000000000000L,0x4000000000000000L,0x2000000000000000L,0x1000000000000000L, 576460752303423488L, 288230376151711744L, 144115188075855872L,  72057594037927936L,  36028797018963968L,  18014398509481984L,   9007199254740992L,   4503599627370496L,   2251799813685248L,   1125899906842624L,    562949953421312L,    281474976710656L,    140737488355328L,     70368744177664L,     35184372088832L,     17592186044416L,      8796093022208L,      4398046511104L,      2199023255552L,      1099511627776L,       549755813888L,       274877906944L,       137438953472L,        68719476736L,        34359738368L,        17179869184L,         8589934592L,         4294967296L,         2147483648L,         1073741824L,          536870912L,          268435456L,          134217728L,           67108864L,           33554432L,           16777216L,            8388608L,            4194304L,            2097152L,            1048576L,             524288L,             262144L,             131072L,              65536L,              32768L,              16384L,               8192L,               4096L,               2048L,               1024L,                512L,                256L,                128L,                 64L,                 32L,                 16L,8L,4L,2L,1L};
const gBytes MR[64]={9223372036854775808L,13835058055282163712L,16140901064495857664L,17293822569102704640L,17870283321406128128L,18158513697557839872L,18302628885633695744L,18374686479671623680L,18410715276690587648L,18428729675200069632L,18437736874454810624L,18442240474082181120L,18444492273895866368L,18445618173802708992L,18446181123756130304L,18446462598732840960L,18446603336221196288L,18446673704965373952L,18446708889337462784L,18446726481523507200L,18446735277616529408L,18446739675663040512L,18446741874686296064L,18446742974197923840L,18446743523953737728L,18446743798831644672L,18446743936270598144L,18446744004990074880L,18446744039349813248L,18446744056529682432L,18446744065119617024L,18446744069414584320L,18446744071562067968L,18446744072635809792L,18446744073172680704L,18446744073441116160L,18446744073575333888L,18446744073642442752L,18446744073675997184L,18446744073692774400L,18446744073701163008L,18446744073705357312L,18446744073707454464L,18446744073708503040L,18446744073709027328L,18446744073709289472L,18446744073709420544L,18446744073709486080L,18446744073709518848L,18446744073709535232L,18446744073709543424L,18446744073709547520L,18446744073709549568L,18446744073709551616L,18446744073709551616L,18446744073709551616L,18446744073709551616L,18446744073709551616L,18446744073709551616L,18446744073709551616L,18446744073709551616L,18446744073709551616L,18446744073709551616L,18446744073709551616L};
#endif
const _unsigned_char mBitsSetTable[256]={
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};

void bitsmove(gBytes *dest,gPos dpos,gBytes *src,gPos spos,gPos npos);
void bitsmoveOR(gBytes *dest,gPos dpos,gBytes *src,gPos spos,gPos npos);
void bitsmoveAND(gBytes *dest,gPos dpos,gBytes *src,gPos spos,gPos npos);
}

using namespace std;
using namespace geco;

void geco::bitsmove(gBytes *dest,gPos dpos,gBytes *src,gPos spos,gPos npos) {
    gPos sptr=spos / mNBits;
    gPos dptr=dpos / mNBits;

    if (npos < 0) {//if(npos < mNBits){
        /*
          if((src+sptr>=dest+dptr)&&(spos>dpos)){
           for(gRelativePos i=0;i<npos;i++){
            if(mGetBit(src,spos+i)) mSetBit(dest,dpos+i);
            else mResetBit(dest,dpos+i);
           }
          }else{
           for(gRelativePos i=(npos-1);i>=0;i--){
            if(mGetBit(src,spos+i)) mSetBit(dest,dpos+i);
            else mResetBit(dest,dpos+i);
           }
          }
        */
    } else {
        gBytes full=~0;
        gSize ncycles = npos / mNBits;
        gSize nrem = npos % mNBits;
        gBytes slbit=spos % mNBits;
        gBytes srbit=mNBits-slbit;
        gBytes dlbit=dpos % mNBits;
        gBytes drbit=mNBits-dlbit;
        gBytes MStart= (drbit==mNBits)?(0):(full << drbit);
        gBytes MEnd = ~MStart;
        gBytes a,m1;

        if ((src+sptr>=dest+dptr)&&(spos>dpos)) {
            if (slbit==dlbit) {
                if (slbit==0) {
                    memmove(dest+dptr,src+sptr,ncycles*sizeof(gBytes));
                    if (nrem>0) {
                        m1=full>>nrem;
                        dest[dptr+ncycles]=(dest[dptr+ncycles] & m1) | (src[sptr+ncycles] & (~m1));
                    }
                } else {
                    if (dlbit+npos<mNBits) {
                        m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                        dest[dptr]=(dest[dptr] & m1) | (src[sptr] & (~m1));
                    } else {
                        gSize lastbits=(dlbit+npos)%mNBits;
                        ncycles=(npos-drbit-lastbits)/mNBits;
                        dest[dptr]=(dest[dptr] & MStart) | (src[sptr] & MEnd);
                        memmove(dest+dptr+1,src+sptr+1,ncycles * mNBytes);
                        m1=(lastbits==mNBits)?(0):(full >> lastbits);
                        dest[dptr+ncycles+1]=(dest[dptr+ncycles+1] & m1) | (src[sptr+ncycles+1] & (~m1));
                    }
                }
            } else {
                if (slbit==0) {
                    for (gPos i=0;i<ncycles;i++) {
                        dest[dptr+i] = (dest[dptr+i] & MStart) | (src[sptr+i] >> dlbit);
                        dest[dptr+i+1] = (dest[dptr+i+1]&MEnd)|(src[sptr+i] << drbit);
                    }
                    if (nrem>0) {
                        if (nrem>drbit) {
                            dest[dptr+ncycles] = (dest[dptr+ncycles] & MStart) | (src[sptr+ncycles] >> dlbit);
                            m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] = (dest[dptr+ncycles+1] & m1) | ((src[sptr+ncycles] << drbit)&(~m1));
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] = (dest[dptr+ncycles] & m1) | ((src[sptr+ncycles] >> dlbit)&(~m1));
                        }
                    }
                } else if (dlbit==0) {
                    for (gPos i=0;i<ncycles;i++) dest[dptr+i]=(src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles]<<slbit);
                        m1=full >> nrem;
                        dest[dptr+ncycles] = (dest[dptr+ncycles] & m1) | (a & (~m1));
                    }
                } else {
                    for (gRelativePos i=0;i< (gRelativePos)ncycles;i++) {
                        a = (src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                        dest[dptr+i] = (dest[dptr+i] & MStart) | (a >> dlbit);
                        dest[dptr+i+1] = (dest[dptr+i+1]&MEnd)|(a << drbit);
                    }
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles] << slbit);
                        if (nrem>drbit) {
                            dest[dptr+ncycles] = (dest[dptr+ncycles] & MStart) | (a >> dlbit);
                            m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] = (dest[dptr+ncycles+1] & m1) | ((a << drbit) & (~m1));
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] = (dest[dptr+ncycles] & m1) | ((a >> dlbit) & (~m1));
                        }
                    }
                }
            }
        } else {
            if (slbit==dlbit) {
                if (slbit==0) {
                    if (nrem>0) {
                        gBytes m1=full>>nrem;
                        dest[dptr+ncycles]=(dest[dptr+ncycles] & m1) | (src[sptr+ncycles] & (~m1));
                    }
                    memmove(dest+dptr,src+sptr,ncycles*sizeof(gBytes));
                } else {
                    if (dlbit+npos<mNBits) {
                        gBytes m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                        dest[dptr]=(dest[dptr] & m1) | (src[sptr] & (~m1));
                    } else {
                        gSize lastbits=(dlbit+npos)%mNBits;
                        ncycles=(npos-drbit-lastbits)/mNBits;
                        m1=(lastbits==mNBits)?(0):(full >> lastbits);
                        dest[dptr+ncycles+1]=(dest[dptr+ncycles+1] & m1) | (src[sptr+ncycles+1] & (~m1));
                        memmove(dest+dptr+1,src+sptr+1,ncycles * mNBytes);
                        dest[dptr]=(dest[dptr] & MStart) | (src[sptr] & MEnd);
                    }
                }
            } else {
                if (slbit==0) {
                    if (nrem>0) {
                        if (nrem>drbit) {
                            gBytes m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] = (dest[dptr+ncycles+1] & m1) | ((src[sptr+ncycles] << drbit)&(~m1));
                            dest[dptr+ncycles] = (dest[dptr+ncycles] & MStart) | (src[sptr+ncycles] >> dlbit);
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] = (dest[dptr+ncycles] & m1) | ((src[sptr+ncycles] >> dlbit)&(~m1));
                        }
                    }
                    for (gRelativePos i=ncycles-1;i>=0;i--) {
                        dest[dptr+i+1] = (dest[dptr+i+1]&MEnd)|(src[sptr+i] << drbit);
                        dest[dptr+i] = (dest[dptr+i] & MStart) | (src[sptr+i] >> dlbit);
                    }
                } else if (dlbit==0) {
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles]<<slbit);
                        m1=full >> nrem;
                        dest[dptr+ncycles] = (dest[dptr+ncycles] & m1) | (a & (~m1));
                    }
                    for (gRelativePos i=(ncycles-1);i>=0;i--) dest[dptr+i]=(src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                } else {
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles] << slbit);
                        if (nrem>drbit) {
                            m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] = (dest[dptr+ncycles+1] & m1 ) | ((a << drbit) & (~m1));
                            dest[dptr+ncycles] = (dest[dptr+ncycles] & MStart) | (a >> dlbit);
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] = (dest[dptr+ncycles] & m1) | ((a >> dlbit) & (~m1));
                        }
                    }
                    for (gRelativePos i=ncycles-1;i>=0;i--) {
                        a = (src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                        dest[dptr+i+1] = (dest[dptr+i+1]&MEnd)|(a << drbit);
                        dest[dptr+i] = (dest[dptr+i] & MStart) | (a >> dlbit);
                    }
                }
            }
        }
    }
}

void geco::bitsmoveOR(gBytes *dest,gPos dpos,gBytes *src,gPos spos,gPos npos) {
    gPos sptr=spos / mNBits;
    gPos dptr=dpos / mNBits;

    if (npos < 0) {//if(npos < mNBits){
        /*
          if((src+sptr>=dest+dptr)&&(spos>dpos)){
           for(gRelativePos i=0;i<npos;i++){
            if(mGetBit(src,spos+i)) mSetBit(dest,dpos+i);
            else mResetBit(dest,dpos+i);
           }
          }else{
           for(gRelativePos i=(npos-1);i>=0;i--){
            if(mGetBit(src,spos+i)) mSetBit(dest,dpos+i);
            else mResetBit(dest,dpos+i);
           }
          }
        */
    } else {
        gBytes full=~0;
        gSize ncycles = npos / mNBits;
        gSize nrem = npos % mNBits;
        gBytes slbit=spos % mNBits;
        gBytes srbit=mNBits-slbit;
        gBytes dlbit=dpos % mNBits;
        gBytes drbit=mNBits-dlbit;
        gBytes MStart= (drbit==mNBits)?(0):(full << drbit);
        gBytes MEnd = ~MStart;
        gBytes a,m1;

        if ((src+sptr>=dest+dptr)&&(spos>dpos)) {
            if (slbit==dlbit) {
                if (slbit==0) {
                    for (gRelativePos i=0;i<(gRelativePos)ncycles;i++) dest[dptr+i]|=src[sptr+i];
                    if (nrem>0) {
                        m1=full>>nrem;
                        dest[dptr+ncycles]|=(dest[dptr+ncycles] & m1) | (src[sptr+ncycles] & (~m1));
                    }
                } else {
                    if (dlbit+npos<mNBits) {
                        m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                        dest[dptr]|=(dest[dptr] & m1) | (src[sptr] & (~m1));
                    } else {
                        gSize lastbits=(dlbit+npos)%mNBits;
                        ncycles=(npos-drbit-lastbits)/mNBits;
                        dest[dptr] |= (dest[dptr] & MStart) | (src[sptr] & MEnd);
                        for (gRelativePos i=1;i<=(gRelativePos)ncycles;i++) dest[dptr+i] |= src[sptr+i];
                        m1=(lastbits==mNBits)?(0):(full >> lastbits);
                        dest[dptr+ncycles+1] |= (dest[dptr+ncycles+1] & m1) | (src[sptr+ncycles+1] & (~m1));
                    }
                }
            } else {
                if (slbit==0) {
                    for (gPos i=0;i<ncycles;i++) {
                        dest[dptr+i] |= (dest[dptr+i] & MStart) | (src[sptr+i] >> dlbit);
                        dest[dptr+i+1] |= (dest[dptr+i+1]&MEnd)|(src[sptr+i] << drbit);
                    }
                    if (nrem>0) {
                        if (nrem>drbit) {
                            dest[dptr+ncycles] |= (dest[dptr+ncycles] & MStart) | (src[sptr+ncycles] >> dlbit);
                            m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] |= (dest[dptr+ncycles+1] & m1) | ((src[sptr+ncycles] << drbit)&(~m1));
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] |= (dest[dptr+ncycles] & m1) | ((src[sptr+ncycles] >> dlbit)&(~m1));
                        }
                    }
                } else if (dlbit==0) {
                    for (gPos i=0;i<ncycles;i++) dest[dptr+i] |= (src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles]<<slbit);
                        m1=full >> nrem;
                        dest[dptr+ncycles] |= (dest[dptr+ncycles] & m1) | (a & (~m1));
                    }
                } else {
                    for (gRelativePos i=0;i<(gRelativePos)ncycles;i++) {
                        a = (src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                        dest[dptr+i] |= (dest[dptr+i] & MStart) | (a >> dlbit);
                        dest[dptr+i+1] |= (dest[dptr+i+1]&MEnd)|(a << drbit);
                    }
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles] << slbit);
                        if (nrem>drbit) {
                            dest[dptr+ncycles] |= (dest[dptr+ncycles] & MStart) | (a >> dlbit);
                            m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] |= (dest[dptr+ncycles+1] & m1) | ((a << drbit) & (~m1));
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] |= (dest[dptr+ncycles] & m1) | ((a >> dlbit) & (~m1));
                        }
                    }
                }
            }
        } else {
            if (slbit==dlbit) {
                if (slbit==0) {
                    if (nrem>0) {
                        gBytes m1=full>>nrem;
                        dest[dptr+ncycles] |= (dest[dptr+ncycles] & m1) | (src[sptr+ncycles] & (~m1));
                    }
                    for (gRelativePos i=(ncycles-1);i>=0;i--) dest[dptr+i] |= src[sptr+i];
                } else {
                    if (dlbit+npos<mNBits) {
                        gBytes m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                        dest[dptr] |= (dest[dptr] & m1) | (src[sptr] & (~m1));
                    } else {
                        gSize lastbits=(dlbit+npos) % mNBits;
                        ncycles=(npos-drbit-lastbits)/mNBits;
                        m1=(lastbits==mNBits)?(0):(full >> lastbits);
                        dest[dptr+ncycles+1] |= (dest[dptr+ncycles+1] & m1) | (src[sptr+ncycles+1] & (~m1));
                        for (gRelativePos i=ncycles;i>0;i--) dest[dptr+i]|=src[sptr+i];
                        dest[dptr] |= (dest[dptr] & MStart) | (src[sptr] & MEnd);
                    }
                }
            } else {
                if (slbit==0) {
                    if (nrem>0) {
                        if (nrem>drbit) {
                            gBytes m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] |= (dest[dptr+ncycles+1] & m1) | ((src[sptr+ncycles] << drbit)&(~m1));
                            dest[dptr+ncycles] |= (dest[dptr+ncycles] & MStart) | (src[sptr+ncycles] >> dlbit);
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] |= (dest[dptr+ncycles] & m1) | ((src[sptr+ncycles] >> dlbit)&(~m1));
                        }
                    }
                    for (gRelativePos i=ncycles-1;i>=0;i--) {
                        dest[dptr+i+1] |= (dest[dptr+i+1]&MEnd)|(src[sptr+i] << drbit);
                        dest[dptr+i] |= (dest[dptr+i] & MStart) | (src[sptr+i] >> dlbit);
                    }
                } else if (dlbit==0) {
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles]<<slbit);
                        m1=full >> nrem;
                        dest[dptr+ncycles] |= (dest[dptr+ncycles] & m1) | (a & (~m1));
                    }
                    for (gRelativePos i=(ncycles-1);i>=0;i--) dest[dptr+i] |= (src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                } else {
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles] << slbit);
                        if (nrem>drbit) {
                            m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] |= (dest[dptr+ncycles+1] & m1 ) | ((a << drbit) & (~m1));
                            dest[dptr+ncycles] |= (dest[dptr+ncycles] & MStart) | (a >> dlbit);
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] |= (dest[dptr+ncycles] & m1) | ((a >> dlbit) & (~m1));
                        }
                    }
                    for (gRelativePos i=ncycles-1;i>=0;i--) {
                        a = (src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                        dest[dptr+i+1] |= (dest[dptr+i+1]&MEnd)|(a << drbit);
                        dest[dptr+i] |= (dest[dptr+i] & MStart) | (a >> dlbit);
                    }
                }
            }
        }
    }
}

void geco::bitsmoveAND(gBytes *dest,gPos dpos,gBytes *src,gPos spos,gPos npos) {
    gPos sptr=spos / mNBits;
    gPos dptr=dpos / mNBits;

    if (npos < 0) {//if(npos < mNBits){
        /*
          if((src+sptr>=dest+dptr)&&(spos>dpos)){
           for(gRelativePos i=0;i<npos;i++){
            if(mGetBit(src,spos+i)) mSetBit(dest,dpos+i);
            else mResetBit(dest,dpos+i);
           }
          }else{
           for(gRelativePos i=(npos-1);i>=0;i--){
            if(mGetBit(src,spos+i)) mSetBit(dest,dpos+i);
            else mResetBit(dest,dpos+i);
           }
          }
        */
    } else {
        gBytes full=~0;
        gSize ncycles = npos / mNBits;
        gSize nrem = npos % mNBits;
        gBytes slbit=spos % mNBits;
        gBytes srbit=mNBits-slbit;
        gBytes dlbit=dpos % mNBits;
        gBytes drbit=mNBits-dlbit;
        gBytes MStart= (drbit==mNBits)?(0):(full << drbit);
        gBytes MEnd = ~MStart;
        gBytes a,m1;

        if ((src+sptr>=dest+dptr)&&(spos>dpos)) {
            if (slbit==dlbit) {
                if (slbit==0) {
                    for (gRelativePos i=0;i<(gRelativePos)ncycles;i++) dest[dptr+i]&=src[sptr+i];
                    if (nrem>0) {
                        m1=full>>nrem;
                        dest[dptr+ncycles] &= (dest[dptr+ncycles] & m1) | (src[sptr+ncycles] & (~m1));
                    }
                } else {
                    if (dlbit+npos<mNBits) {
                        m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                        dest[dptr] &= (dest[dptr] & m1) | (src[sptr] & (~m1));
                    } else {
                        gSize lastbits=(dlbit+npos)%mNBits;
                        ncycles=(npos-drbit-lastbits)/mNBits;
                        dest[dptr] &= (dest[dptr] & MStart) | (src[sptr] & MEnd);
                        for (gRelativePos i=1;i<=(gRelativePos)ncycles;i++) dest[dptr+i] &= src[sptr+i];
                        m1=(lastbits==mNBits)?(0):(full >> lastbits);
                        dest[dptr+ncycles+1] &= (dest[dptr+ncycles+1] & m1) | (src[sptr+ncycles+1] & (~m1));
                    }
                }
            } else {
                if (slbit==0) {
                    for (gPos i=0;i<ncycles;i++) {
                        dest[dptr+i] &= (dest[dptr+i] & MStart) | (src[sptr+i] >> dlbit);
                        dest[dptr+i+1] &= (dest[dptr+i+1]&MEnd)|(src[sptr+i] << drbit);
                    }
                    if (nrem>0) {
                        if (nrem>drbit) {
                            dest[dptr+ncycles] &= (dest[dptr+ncycles] & MStart) | (src[sptr+ncycles] >> dlbit);
                            m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] &= (dest[dptr+ncycles+1] & m1) | ((src[sptr+ncycles] << drbit)&(~m1));
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] &= (dest[dptr+ncycles] & m1) | ((src[sptr+ncycles] >> dlbit)&(~m1));
                        }
                    }
                } else if (dlbit==0) {
                    for (gPos i=0;i<ncycles;i++) dest[dptr+i] &= (src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles]<<slbit);
                        m1=full >> nrem;
                        dest[dptr+ncycles] &= (dest[dptr+ncycles] & m1) | (a & (~m1));
                    }
                } else {
                    for (gRelativePos i=0;i<(gRelativePos)ncycles;i++) {
                        a = (src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                        dest[dptr+i] &= (dest[dptr+i] & MStart) | (a >> dlbit);
                        dest[dptr+i+1] &= (dest[dptr+i+1]&MEnd)|(a << drbit);
                    }
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles] << slbit);
                        if (nrem>drbit) {
                            dest[dptr+ncycles] &= (dest[dptr+ncycles] & MStart) | (a >> dlbit);
                            m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] &= (dest[dptr+ncycles+1] & m1) | ((a << drbit) & (~m1));
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] &= (dest[dptr+ncycles] & m1) | ((a >> dlbit) & (~m1));
                        }
                    }
                }
            }
        } else {
            if (slbit==dlbit) {
                if (slbit==0) {
                    if (nrem>0) {
                        gBytes m1=full>>nrem;
                        dest[dptr+ncycles] &= (dest[dptr+ncycles] & m1) | (src[sptr+ncycles] & (~m1));
                    }
                    for (gRelativePos i=(ncycles-1);i>=0;i--) dest[dptr+i] &= src[sptr+i];
                } else {
                    if (dlbit+npos<mNBits) {
                        gBytes m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                        dest[dptr] &= (dest[dptr] & m1) | (src[sptr] & (~m1));
                    } else {
                        gSize lastbits=(dlbit+npos) % mNBits;
                        ncycles=(npos-drbit-lastbits)/mNBits;
                        m1=(lastbits==mNBits)?(0):(full >> lastbits);
                        dest[dptr+ncycles+1] &= (dest[dptr+ncycles+1] & m1) | (src[sptr+ncycles+1] & (~m1));
                        for (gRelativePos i=ncycles;i>0;i--) dest[dptr+i] &= src[sptr+i];
                        dest[dptr] &= (dest[dptr] & MStart) | (src[sptr] & MEnd);
                    }
                }
            } else {
                if (slbit==0) {
                    if (nrem>0) {
                        if (nrem>drbit) {
                            gBytes m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] &= (dest[dptr+ncycles+1] & m1) | ((src[sptr+ncycles] << drbit)&(~m1));
                            dest[dptr+ncycles] &= (dest[dptr+ncycles] & MStart) | (src[sptr+ncycles] >> dlbit);
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] &= (dest[dptr+ncycles] & m1) | ((src[sptr+ncycles] >> dlbit)&(~m1));
                        }
                    }
                    for (gRelativePos i=ncycles-1;i>=0;i--) {
                        dest[dptr+i+1] &= (dest[dptr+i+1]&MEnd)|(src[sptr+i] << drbit);
                        dest[dptr+i] &= (dest[dptr+i] & MStart) | (src[sptr+i] >> dlbit);
                    }
                } else if (dlbit==0) {
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles]<<slbit);
                        m1=full >> nrem;
                        dest[dptr+ncycles] &= (dest[dptr+ncycles] & m1) | (a & (~m1));
                    }
                    for (gRelativePos i=(ncycles-1);i>=0;i--) dest[dptr+i] &= (src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                } else {
                    if (nrem>0) {
                        if (nrem>srbit) a=(src[sptr+ncycles] << slbit) | (src[sptr+ncycles+1] >> srbit);
                        else a=(src[sptr+ncycles] << slbit);
                        if (nrem>drbit) {
                            m1= full >> (nrem-drbit);
                            dest[dptr+ncycles+1] &= (dest[dptr+ncycles+1] & m1 ) | ((a << drbit) & (~m1));
                            dest[dptr+ncycles] &= (dest[dptr+ncycles] & MStart) | (a >> dlbit);
                        } else {
                            m1= ((drbit==mNBits)?(0):(full << drbit)) | ((dlbit+nrem==mNBits)?(0):(full >> (dlbit+nrem)));
                            dest[dptr+ncycles] &= (dest[dptr+ncycles] & m1) | ((a >> dlbit) & (~m1));
                        }
                    }
                    for (gRelativePos i=ncycles-1;i>=0;i--) {
                        a = (src[sptr+i] << slbit) | (src[sptr+i+1] >> srbit);
                        dest[dptr+i+1] &= (dest[dptr+i+1]&MEnd)|(a << drbit);
                        dest[dptr+i] &= (dest[dptr+i] & MStart) | (a >> dlbit);
                    }
                }
            }
        }
    }
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// class gBitsArray definition
//--------------------------------------------------------------------------
/** @brief Empty Constructor
 *
 * Instatiates an empty gBitsArray object
 */
gBitsArray::gBitsArray() {
    i_val=0;
    i_alength=0;
    i_mlength=0;
}

/** @brief Non initializing length constructor
 *
 * This constructor builds a bits array of the given lengrth
 * without initialization.
 * @param length The desired array length
 */
gBitsArray::gBitsArray(gSize length) {
    if (length>0) {
        i_alength=length;
        i_mlength=length/mNBits+1;
        i_val=new gBytes[i_mlength];
        i_val[0]=0;
        i_val[i_mlength-1]=0;
    } else {
        i_val=0;
        i_alength=0;
        i_mlength=0;
    }
}

/** @brief Initializing length constructor.
 *
 * This constructor builds a bits array of the given length
 * initializing all its element to the boolean value specified.
 * @param length The desired array length
 * @param value The initialization value
 */
gBitsArray::gBitsArray(gSize length,bool value) {
    if (length>0) {
        i_alength=length;
        i_mlength=length/mNBits+1;
        i_val=new gBytes[i_mlength];
        memset(i_val,(value)?(255):(0),mNBytes*i_mlength);
    } else {
        i_val=0;
        i_alength=0;
        i_mlength=0;
    }
}

/** @brief Range copy constructor
 *
 * This is a range copy contructor, it builds a new object using values from an existing one in the range \a start - \a end
 * @param array The array to copy from
 * @param start The zero-based position from which values have to be copied. Defaults to 0 (beginning).
 * @param end The one-based position to which values have to be copied. A value of 0 (default) makes the array to be copied to the end.
 */
gBitsArray::gBitsArray(const gBitsArray & array,gPos start, gPos end) {
    gPos aend=(end==0)?(array.i_alength):(end);
    if ((start==0)&&(aend==array.i_alength)) {
        i_val = new gBytes[array.i_mlength];
        i_mlength=array.i_mlength;
        i_alength=array.i_alength;
        memmove(i_val,array.i_val,array.i_mlength*mNBytes);
    } else {
        i_val=0;
        i_alength=0;
        i_mlength=0;
        gPos aend=(end==0)?(array.i_alength):(end);
        copyValues(0,array,start,aend);
    }

}

/** @brief Destructor */
gBitsArray::~gBitsArray() {
    if (i_val) delete [] i_val;
    i_val=0;
    i_alength=0;
    i_mlength=0;
}

/** @brief Bits array length
 *
 * Return the size in bits of the current array
 * @return the size of this array
 */
gSize gBitsArray::getSize() const {
    return i_alength;
}

/** @brief Set bits counter
 *
 * This method returns the number of bits set in the specified range.
 * @param start Range start (zero based)
 * @param end Range end (one based)
 * @return number of set bits
 */
gSize gBitsArray::getSetBitsCount(gPos start,gPos end) const {
    gSize count=0,i;
    for(gPos i=start;i<end;i++) count+=(gSize) getValue(i);
    /* This following is the smart way baut it doesn't work!!!!!
    unsigned char *cptr;
    gBytes val;
    gRelativePos a = start / mNBits;
    gRelativePos b = (end) / mNBits;
    if (a==b) {
        val=i_val[a] & (((~0) >> (start % mNBits)) << ( mNBits - (end % mNBits)));
        cptr = ((unsigned char *) &val);
        for (i=0;i<mNBytes;i++) {
            if (*cptr!=0) count+=mBitsSetTable[*cptr];
            cptr++;
        }
    } else {
        val=(i_val[a] & ((~0) >> (start % mNBits)));
        if (val>0) {
            cptr = ((unsigned char *) & val);
            for (i=0;i<mNBytes;i++) {
                count+=mBitsSetTable[*cptr];
                cptr++;
            }
        }
        cptr = (unsigned char *)(i_val+a+1);
        gSize num=(b-a-1)*mNBytes;
        for (i=0;i<num;i++) {
            count+=mBitsSetTable[*cptr];
            cptr++;
        }
        val=i_val[b] & ((~0) << ( mNBits - (end % mNBits)));
        if (val>0) {
            cptr = ((unsigned char *) &val);
            for (i=0;i<mNBytes;i++) {
                count+=mBitsSetTable[*cptr];
                cptr++;
            }
        }
    }
    */ 
    return count;
}

/** @brief Change array size
 *
 * this method changes the array size to the desired \a size optionally initializing added positions
 * to a specified value
 * @param size The desired array size
 * @param initvalues Whether or not to initialize added positions to the value specified by fillvalue. Defaults to true.
 * @param fillvalue The value to which initialize added positions (ignored if initvalues is false). Defaults to false.
 */
void gBitsArray::setSize(gSize size,bool initvalues,bool fillvalue) {
    if (size==0) {
        delete [] i_val;
        i_val=0;
        i_mlength=0;
        i_alength=0;
    } else {
        gSize mlen=size/mNBits+1;
        if (size!=i_alength) {
            gBytes *nmem=new gBytes[mlen];
            if (i_alength==0) {
                if (initvalues) memset(nmem,(fillvalue)?(255):(0),mNBytes*mlen);
            } else {
                if (size>i_alength) {
                    //if(initvalues) memset(nmem+i_mlength-1,(fillvalue)?(255):(0),mNBytes*(mlen-i_mlength));
                    if (initvalues) memset(nmem,(fillvalue)?(255):(0),mNBytes*mlen);
                    bitsmove(nmem,0,i_val,0,i_alength);
                } else {
                    if (initvalues) memset(nmem,(fillvalue)?(255):(0),mNBytes*mlen);
                    bitsmove(nmem,0,i_val,0,size);
                }
                delete [] i_val;
            }
            i_val=nmem;
            i_mlength=mlen;
            i_alength=size;
        }
    }
}

/** @brief  Copy bit values from a range
 *
 * This method copies aend-astart bits from another bits array into this one starting from start
 * @param start Destination start position (zero based)
 * @param array Source array
 * @param astart Start position in the source array (zero based)
 * @param aend End position in the source array (one based)
 */
void gBitsArray::copyValues(gPos start,const gBitsArray & array,gPos astart,gPos aend) {
    if (start+aend-astart > i_alength) setSize(start+(aend-astart),false);
    bitsmove(i_val,start,array.i_val,astart,aend-astart);
}

/** @brief Ranged or with another array
 *
 * Set aend-astart bits starting from start to the values given by or-ing them with
 * corresponding bits in the providend array between astart and aend.
 * @param start Destination starting position
 * @param array Suorce array
 * @param astart Starting position in the provided array. (zero based)
 * @param aend Ending position in the provided array (one based)
 */
void gBitsArray::Or(gPos start,const gBitsArray & array,gPos astart,gPos aend) {
    if (start+aend-astart > i_alength) setSize(start+(aend-astart),false);
    bitsmoveOR(i_val,start,array.i_val,astart,aend-astart);
}

/** @brief Ranged and with another array
 *
 * Set aend-astart bits starting from start to the values given by and-ing them with
 * corresponding bits in the providend array between astart and aend.
 * @param start destination starting position
 * @param array provided bits array
 * @param astart Starting position in the provided array (zero based)
 * @param aend Ending pèosition in the provided array (one based)
 */
void gBitsArray::And(gPos start,const gBitsArray & array,gPos astart,gPos aend) {
    if (start+aend-astart > i_alength) setSize(start+(aend-astart),false);
    bitsmoveAND(i_val,start,array.i_val,astart,aend-astart);
}

/** @brief Set all bits
 */
void gBitsArray::setAll() {
    memset(i_val,255,mNBytes*i_mlength);
}

/** @brief Reset all bits
 */
void gBitsArray::resetAll() {
    memset(i_val,0,mNBytes*i_mlength);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// class gException
//--------------------------------------------------------------------------
/** @brief Reason contsructor
 *
 * @param reason: string containing the event description.
 */
gException::gException(const char *reason):exception() {
    if (reason) i_reason=reason;
}

/** @brief Copy contructor
 *
 * @param e: the object to copy from
 */
gException::gException(const gException & e):exception() {
    *this=e;
}

/** @brief Throwable destructor */
gException::~gException() throw() {
}

/** @brief Assignment operator
 *
 * @param e: the object to copy from
 * @return a reference to this
 */
gException & gException::operator = (const gException &e) {
    i_reason=e.i_reason;
    return *this;
}

/** @brief Get event description
 *
 * This method return a string containing
 * a description of the event that caused the
 * exception to be thrown.
 * @return a const char * containing the event description.
 */
const char* gException::what() const throw() {
    return i_reason.c_str();
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// class gRetrieverImplementation definitions
//--------------------------------------------------------------------------
/** @brief Default constructor
 * 
 * This constructor initialize the implementation refCount mechanism: it must be always
 * called from derived objects
*/
gRetrieverImplementation::gRetrieverImplementation(){
  i_refCount=0;
}

/** @brief Virtual destructor
 * 
 * Must be reimplemented by any derived object
*/
gRetrieverImplementation::~gRetrieverImplementation(){
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// class gRetriever definitions
//--------------------------------------------------------------------------
/** @brief Default constructor
 * 
 * This constructor initialize the refCount mechanism: it must be always
 * called from derived objects
*/
gRetriever::gRetriever(){
  i_implementation=NULL;
}

/** @brief Implementation constructor
 * 
 * This constructor initializes the refCount mechanism: it must be always
 * called from derived objects
 * @param implementation provide the retriever implementation
*/
gRetriever::gRetriever( const gRetrieverImplementation & implementation ){
  i_implementation=implementation.clone();
  i_implementation->i_refCount++;
}

/** @brief Copy constructor
 * 
 * This constructor initializes the refCount mechanism: it must be always
 * called from derived objects
 * @param retriever provide the retriever to copy from
*/
gRetriever::gRetriever( const gRetriever & retriever){
  i_implementation=retriever.i_implementation;
  i_implementation->i_refCount++;
}

/** @brief Desctructor */
gRetriever::~gRetriever(){
  if(i_implementation){
    i_implementation->i_refCount--;
    if(i_implementation->i_refCount==0) delete i_implementation;
  }
}

gRetriever & gRetriever::operator = (const gRetriever & retriever){
  if(i_implementation){
    i_implementation->i_refCount--;
    if(i_implementation->i_refCount==0) delete i_implementation;
  }
  i_implementation=retriever.i_implementation;
  i_implementation->i_refCount++;
}
/** @brief Implementation access 
 *
 * This protected member function allow derived classes to access the 
 * implementation of the retreiver. Usually they should cast this 
 * reference to the sepcific type in  order to acces its methods
 * @return A const gRetrieverImplementation reference
 */
const gRetrieverImplementation & gRetriever::getImplementation() const { 
  return * i_implementation;
}
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
// template class gArrayInternal definition
//--------------------------------------------------------------------------
// template class em_bas_array definition is in the gArray.cpp file in order
// to be available to the gArray template explicit instantiation
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void geco_running() {
}
//--------------------------------------------------------------------------


