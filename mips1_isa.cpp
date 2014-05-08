/**
 * @file      mips1_isa.cpp
 * @author    Sandro Rigo
 *            Marcus Bartholomeu
 *            Alexandro Baldassin (acasm information)
 *
 *            The ArchC Team
 *            http://www.archc.org/
 *
 *            Computer Systems Laboratory (LSC)
 *            IC-UNICAMP
 *            http://www.lsc.ic.unicamp.br/
 *
 * @version   1.0
 * @date      Mon, 19 Jun 2006 15:50:52 -0300
 * 
 * @brief     The ArchC i8051 functional model.
 * 
 * @attention Copyright (C) 2002-2006 --- The ArchC Team
 * 
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version. 
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 * 
 * You should have received a copy of the GNU General Public License 
 * along with this program; if not, write to the Free Software 
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include  "mips1_isa.H"
#include  "mips1_isa_init.cpp"
#include  "mips1_bhv_macros.H"
#include  <map>
#include  <vector>

//If you want debug information for this model, uncomment next line
#define DEBUG_MODEL
#include "ac_debug_model.H"


//!User defined macros to reference registers.
#define Ra 31
#define Sp 29

// 'using namespace' statement to allow access to all
// mips1-specific datatypes
using namespace mips1_parms;

//counters 5-Stage
int countForward = 0;
int countStalls = 0;
int branchStalls = 0;

//counters 7-Stage
int countForward7 = 0;
int countStalls7 = 0;
int branchStalls7 = 0;

//7-Stage Pipeline controls
std::vector <std::vector <int> > stage7 (3,std::vector<int>(7,0));
std::vector <std::vector <bool> > stage7ctrl (2,std::vector<bool>(7,0));

//5-Stage Pipeline controls
int rdIfId = 0, rdIdEx = 0, rdExMem = 0, rdMemWb = 0;
int rtIfId = 0, rtIdEx = 0, rtExMem = 0, rtMemWb = 0;
int rsIfId = 0, rsIdEx = 0, rsExMem = 0, rsMemWb = 0;
bool preMemRead = 0, Id_Ex_MemRead = 0, Ex_Mem_MemRead = 0;
bool preRegWrite = 0, Id_Ex_RegWrite = 0, Ex_Mem_RegWrite = 0, Mem_Wb_RegWrite = 0;

//Dynamic Branch Prediction
std::map<int,int> bp;
std::map<int,int>::iterator it;
int dynamicStalls = 0;
int dynamicStalls7 = 0;
std::map<int,std::vector<int> > bhist;
std::map<int,std::vector<int> >::iterator ithist;
int dynamicStalls2 = 0;
int dynamicStalls27 = 0;
int predict;
void movepipe()
{
  Mem_Wb_RegWrite = Ex_Mem_RegWrite;
  Ex_Mem_RegWrite = Id_Ex_RegWrite;
  Id_Ex_RegWrite = preRegWrite;
  
  Ex_Mem_MemRead = Id_Ex_MemRead;
  Id_Ex_MemRead = preMemRead;
  
  rdMemWb = rdExMem;
  rdExMem = rdIdEx;
  rdIdEx = rdIfId;
  
  rtMemWb = rtExMem;
  rtExMem = rtIdEx;
  rtIdEx = rtIfId;
  
  rsMemWb = rsExMem;
  rsExMem = rsIdEx;
  rsIdEx = rsIfId;

  stage7[0].pop_back();//Rd
  stage7[1].pop_back();//Rt
  stage7[2].pop_back();//Rs

  stage7[0].insert(stage7[0].begin(), 0);//Rd
  stage7[1].insert(stage7[1].begin(), 0);//Rt
  stage7[2].insert(stage7[2].begin(), 0);//Rs

  stage7ctrl[0].pop_back();//RegWrite
  stage7ctrl[1].pop_back();//MemRead

  stage7ctrl[0].insert(stage7ctrl[0].begin(), 0);//RegWrite
  stage7ctrl[1].insert(stage7ctrl[1].begin(), 0);//MemRead
}


void forward_check(){
  int flag=0;
  int flag7=0;
  //Rs
  if (Ex_Mem_RegWrite && rdExMem != 0 && rdExMem == rsIdEx){//Ex hazard
    countForward+=2;
    flag=2;
  } else if (Mem_Wb_RegWrite && rdMemWb != 0 && rdMemWb == rsIdEx) {//Mem hazard
    countForward++;
    flag=1;
  }
  //Rt
  if (Ex_Mem_RegWrite && rdExMem != 0 && rdExMem == rtIdEx){//Ex hazard
    if (flag == 0) countForward+=2; else if (flag == 1) countForward++;
  } else if (Mem_Wb_RegWrite && rdMemWb != 0 && rdMemWb == rtIdEx) {//Mem hazard
    if (flag == 0) countForward++;
  }
  
  //Rs7
  if (stage7ctrl[0][3] && stage7[0][3] != 0 && stage7[0][3] == stage7[2][2]){//Ex hazard
    countForward7+=3;
    flag7=3;
  } else if (stage7ctrl[0][4] && stage7[0][4] != 0 && stage7[0][4] == stage7[2][2]){//Mt hazard
    countForward7+=2;
    flag7=2;
  } else if (stage7ctrl[0][5] && stage7[0][5] != 0 && stage7[0][5] == stage7[2][2]){//Mem hazard
    countForward7+=1;
    flag7=1;
  }
  //Rt7
  if (stage7ctrl[0][3] && stage7[0][3] != 0 && stage7[0][3] == stage7[1][2]){//Ex hazard
    if (flag7 == 0) countForward7+=3; else if (flag7 == 1) countForward7+=2; else if (flag7 == 2) countForward7+=1;
  } else if (stage7ctrl[0][4] && stage7[0][4] != 0 && stage7[0][4] == stage7[1][2]){//Mt hazard
    if (flag7 == 0) countForward7+=2; else if (flag7 == 1) countForward7+=1;
  } else if (stage7ctrl[0][5] && stage7[0][5] != 0 && stage7[0][5] == stage7[1][2]){//Mem hazard
    if (flag7 == 0) countForward7+=1;
  }
}

void stall_check(){
  
  if (Id_Ex_MemRead && (rtIdEx == rsIfId || rtIdEx == rtIfId))
    countStalls++;
  
  if (stage7ctrl[1][2] && (stage7[1][2] == stage7[2][1] || stage7[1][2] == stage7[1][1]))
    countStalls7+=2;
  if (stage7ctrl[1][2] && (stage7[1][2] == stage7[2][0] || stage7[1][2] == stage7[1][0]))
    countStalls7++;

}

void pipecontrol5(){
  int i;
  
  for (i = 0; i < 6; i++){
    movepipe();
    preRegWrite = 0;
    preMemRead = 0;
    rdIfId = 0;
    rtIfId = 0;
    rsIfId = 0;
    stall_check();
    forward_check();
  }
  
  //dbg_printf("rd: %2d %2d %2d %2d\nrt: %2d %2d %2d %2d\nrs: %2d %2d %2d %2d\n", rdIfId, rdIdEx, rdExMem, rdMemWb
  //, rtIfId, rtIdEx, rtExMem, rtMemWb, rsIfId, rsIdEx, rsExMem, rsMemWb);
  
}

//!Generic instruction behavior method.
void ac_behavior( instruction )
{ 
  dbg_printf("----- PC=%#x ----- %lld\n", (int) ac_pc, ac_instr_counter);
  //  dbg_printf("----- PC=%#x NPC=%#x ----- %lld\n", (int) ac_pc, (int)npc, ac_instr_counter);
#ifndef NO_NEED_PC_UPDATE
  ac_pc = npc;
  npc = ac_pc + 4;
#endif 
  printf("2 %x din\n", (int)ac_pc);
  movepipe();
};
 
//! Instruction Format behavior methods.
void ac_behavior( Type_R ){
  //5stage
  preRegWrite = 1;
  preMemRead = 0;
  rdIfId = rd;
  rtIfId = rt;
  rsIfId = rs;
  
  //7stage
  stage7ctrl[0][0] = 1;
  stage7ctrl[1][0] = 0;
  stage7[0][0] = rd;
  stage7[1][0] = rt;
  stage7[2][0] = rs;
  
  stall_check();
  forward_check();
  //5stage_end
}
void ac_behavior( Type_R_Jump ){
  //5stage
  preRegWrite = 1;
  preMemRead = 0;
  rdIfId = rd;
  rtIfId = rt;
  rsIfId = rs;
  
  //7stage
  stage7ctrl[0][0] = 1;
  stage7ctrl[1][0] = 0;
  stage7[0][0] = rd;
  stage7[1][0] = rt;
  stage7[2][0] = rs;
  
  stall_check();
  forward_check();
  //5stage_end
}
void ac_behavior( Type_NOP ){
  //5stage
  preRegWrite = 0;
  preMemRead = 0;
  rdIfId = 0;
  rtIfId = 0;
  rsIfId = 0;
  
  //7stage
  stage7ctrl[0][0] = 0;
  stage7ctrl[1][0] = 0;
  stage7[0][0] = 0;
  stage7[1][0] = 0;
  stage7[2][0] = 0;
  
  stall_check();
  forward_check();
  //5stage_end
}
void ac_behavior( Type_I ){
  //5stage
  preRegWrite = 0;
  preMemRead = 0;
  rdIfId = 0;
  rtIfId = rt;
  rsIfId = rs;
  
  //7stage
  stage7ctrl[0][0] = 0;
  stage7ctrl[1][0] = 0;
  stage7[0][0] = 0;
  stage7[1][0] = rt;
  stage7[2][0] = rs;
  
  stall_check();
  forward_check();
  //5stage_end
  
  //Verify if branch exist in hash
  it = bp.find((int)ac_pc);
  if(it == bp.end()){
    bp.insert(std::pair<int,int>((int) ac_pc, 0));
  }

  ithist = bhist.find((int)ac_pc);
  if(ithist == bhist.end()){
      bhist.insert(std::pair<int, std::vector<int> >((int)ac_pc, std::vector<int> (4, 0)));
  } else {
      predict = bhist[(int)ac_pc][bp[(int)ac_pc]];
  }
}
void ac_behavior( Type_I_STORE ){
  //5stage
  preRegWrite = 0;
  preMemRead = 0;
  rdIfId = 0;
  rtIfId = rt;
  rsIfId = rs;
  
  //7stage
  stage7ctrl[0][0] = 0;
  stage7ctrl[1][0] = 0;
  stage7[0][0] = 0;
  stage7[1][0] = rt;
  stage7[2][0] = rs;
  
  printf("1 %x din\n", RB[rs]+imm);
  
  stall_check();
  forward_check();
  //5stage_end
}
void ac_behavior( Type_I_LOAD ){
  //5stage
  preRegWrite = 1;
  preMemRead = 0;
  rdIfId = rt;
  rtIfId = rt;
  rsIfId = rs;
  
  //7stage
  stage7ctrl[0][0] = 1;
  stage7ctrl[1][0] = 0;
  stage7[0][0] = rt;
  stage7[1][0] = rt;
  stage7[2][0] = rs;
  
  stall_check();
  forward_check();
  //5stage_end
}
void ac_behavior( Type_I_LOAD_MEM ){
  //5stage
  preRegWrite = 1;
  preMemRead = 1;
  rdIfId = rt;
  rtIfId = rt;
  rsIfId = rs;
  
  //7stage
  stage7ctrl[0][0] = 1;
  stage7ctrl[1][0] = 1;
  stage7[0][0] = rt;
  stage7[1][0] = rt;
  stage7[2][0] = rs;
  
  printf("0 %x din\n", RB[rs]+imm);
  
  stall_check();
  forward_check();
  //5stage_end
}
void ac_behavior( Type_J ){
  //5stage
  preRegWrite = 0;
  preMemRead = 0;
  rdIfId = 0;
  rtIfId = 0;
  rsIfId = 0;
  
  //7stage
  stage7ctrl[0][0] = 0;
  stage7ctrl[1][0] = 0;
  stage7[0][0] = 0;
  stage7[1][0] = 0;
  stage7[2][0] = 0;
  
  stall_check();
  forward_check();
  //5stage_end
}
 
//!Behavior called before starting simulation
void ac_behavior(begin)
{
  dbg_printf("@@@ begin behavior @@@\n");
  RB[0] = 0;
  npc = ac_pc + 4;

  // Is is not required by the architecture, but makes debug really easier
  for (int regNum = 0; regNum < 32; regNum ++)
    RB[regNum] = 0;
  hi = 0;
  lo = 0;
}

//!Behavior called after finishing simulation
void ac_behavior(end)
{
  pipecontrol5();
  printf("--- Pipeline 5 ---\n");
  printf("$$$ Stalls: %d $$$\n", countStalls);
  printf("$$$ Forwards: %d $$$\n", countForward);
  printf("$$$ BranchStalls: %d $$$\n", branchStalls);

  printf("--- Pipeline 7 ---\n");
  printf("$$$ Stalls: %d $$$\n", countStalls7);
  printf("$$$ Forwards: %d $$$\n", countForward7);
  printf("$$$ BranchStalls: %d $$$\n", branchStalls7);

  printf("--- Branch Prediction ---\n");
  printf("$$$ DynamicStalls1: %d $$$\n", dynamicStalls);
  printf("$$$ DynamicStalls1 - 7: %d $$$\n", dynamicStalls7);
  printf("$$$ DynamicStalls2: %d $$$\n", dynamicStalls2);
  printf("$$$ DynamicStalls2 - 7: %d $$$\n", dynamicStalls27);
  
  dbg_printf("@@@ end behavior @@@\n");
}


//!Instruction lb behavior method.
void ac_behavior( lb )
{
  char byte;
  dbg_printf("lb r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  byte = DM.read_byte(RB[rs]+ imm);
  RB[rt] = (ac_Sword)byte ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lbu behavior method.
void ac_behavior( lbu )
{
  unsigned char byte;
  dbg_printf("lbu r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  byte = DM.read_byte(RB[rs]+ imm);
  RB[rt] = byte ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lh behavior method.
void ac_behavior( lh )
{
  short int half;
  dbg_printf("lh r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  half = DM.read_half(RB[rs]+ imm);
  RB[rt] = (ac_Sword)half ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lhu behavior method.
void ac_behavior( lhu )
{
  unsigned short int  half;
  half = DM.read_half(RB[rs]+ imm);
  RB[rt] = half ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lw behavior method.
void ac_behavior( lw )
{
  dbg_printf("lw r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  RB[rt] = DM.read(RB[rs]+ imm);
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lwl behavior method.
void ac_behavior( lwl )
{
  dbg_printf("lwl r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  unsigned int addr, offset;
  ac_Uword data;

  addr = RB[rs] + imm;
  offset = (addr & 0x3) * 8;
  data = DM.read(addr & 0xFFFFFFFC);
  data <<= offset;
  data |= RB[rt] & ((1<<offset)-1);
  RB[rt] = data;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lwr behavior method.
void ac_behavior( lwr )
{
  dbg_printf("lwr r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  unsigned int addr, offset;
  ac_Uword data;

  addr = RB[rs] + imm;
  offset = (3 - (addr & 0x3)) * 8;
  data = DM.read(addr & 0xFFFFFFFC);
  data >>= offset;
  data |= RB[rt] & (0xFFFFFFFF << (32-offset));
  RB[rt] = data;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction sb behavior method.
void ac_behavior( sb )
{
  unsigned char byte;
  dbg_printf("sb r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  byte = RB[rt] & 0xFF;
  DM.write_byte(RB[rs] + imm, byte);
  dbg_printf("Result = %#x\n", (int) byte);
};

//!Instruction sh behavior method.
void ac_behavior( sh )
{
  unsigned short int half;
  dbg_printf("sh r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  half = RB[rt] & 0xFFFF;
  DM.write_half(RB[rs] + imm, half);
  dbg_printf("Result = %#x\n", (int) half);
};

//!Instruction sw behavior method.
void ac_behavior( sw )
{
  dbg_printf("sw r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  DM.write(RB[rs] + imm, RB[rt]);
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction swl behavior method.
void ac_behavior( swl )
{
  dbg_printf("swl r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  unsigned int addr, offset;
  ac_Uword data;

  addr = RB[rs] + imm;
  offset = (addr & 0x3) * 8;
  data = RB[rt];
  data >>= offset;
  data |= DM.read(addr & 0xFFFFFFFC) & (0xFFFFFFFF << (32-offset));
  DM.write(addr & 0xFFFFFFFC, data);
  dbg_printf("Result = %#x\n", data);
};

//!Instruction swr behavior method.
void ac_behavior( swr )
{
  dbg_printf("swr r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  unsigned int addr, offset;
  ac_Uword data;

  addr = RB[rs] + imm;
  offset = (3 - (addr & 0x3)) * 8;
  data = RB[rt];
  data <<= offset;
  data |= DM.read(addr & 0xFFFFFFFC) & ((1<<offset)-1);
  DM.write(addr & 0xFFFFFFFC, data);
  dbg_printf("Result = %#x\n", data);
};

//!Instruction addi behavior method.
void ac_behavior( addi )
{
  dbg_printf("addi r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  RB[rt] = RB[rs] + imm;
  dbg_printf("Result = %#x\n", RB[rt]);
  //Test overflow
  if ( ((RB[rs] & 0x80000000) == (imm & 0x80000000)) &&
       ((imm & 0x80000000) != (RB[rt] & 0x80000000)) ) {
    fprintf(stderr, "EXCEPTION(addi): integer overflow.\n"); exit(EXIT_FAILURE);
  }
};

//!Instruction addiu behavior method.
void ac_behavior( addiu )
{
  dbg_printf("addiu r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  RB[rt] = RB[rs] + imm;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction slti behavior method.
void ac_behavior( slti )
{
  dbg_printf("slti r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  // Set the RD if RS< IMM
  if( (ac_Sword) RB[rs] < (ac_Sword) imm )
    RB[rt] = 1;
  // Else reset RD
  else
    RB[rt] = 0;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction sltiu behavior method.
void ac_behavior( sltiu )
{
  dbg_printf("sltiu r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  // Set the RD if RS< IMM
  if( (ac_Uword) RB[rs] < (ac_Uword) imm )
    RB[rt] = 1;
  // Else reset RD
  else
    RB[rt] = 0;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction andi behavior method.
void ac_behavior( andi )
{	
  dbg_printf("andi r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  RB[rt] = RB[rs] & (imm & 0xFFFF) ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction ori behavior method.
void ac_behavior( ori )
{	
  dbg_printf("ori r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  RB[rt] = RB[rs] | (imm & 0xFFFF) ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction xori behavior method.
void ac_behavior( xori )
{	
  dbg_printf("xori r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  RB[rt] = RB[rs] ^ (imm & 0xFFFF) ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lui behavior method.
void ac_behavior( lui )
{	
  dbg_printf("lui r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  // Load a constant in the upper 16 bits of a register
  // To achieve the desired behaviour, the constant was shifted 16 bits left
  // and moved to the target register ( rt )
  RB[rt] = imm << 16;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction add behavior method.
void ac_behavior( add )
{
  dbg_printf("add r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] + RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
  //Test overflow
  if ( ((RB[rs] & 0x80000000) == (RB[rd] & 0x80000000)) &&
       ((RB[rd] & 0x80000000) != (RB[rt] & 0x80000000)) ) {
    fprintf(stderr, "EXCEPTION(add): integer overflow.\n"); exit(EXIT_FAILURE);
  }
};

//!Instruction addu behavior method.
void ac_behavior( addu )
{
  dbg_printf("addu r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] + RB[rt];
  //cout << "  RS: " << (unsigned int)RB[rs] << " RT: " << (unsigned int)RB[rt] << endl;
  //cout << "  Result =  " <<  (unsigned int)RB[rd] <<endl;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction sub behavior method.
void ac_behavior( sub )
{
  dbg_printf("sub r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] - RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
  //TODO: test integer overflow exception for sub
};

//!Instruction subu behavior method.
void ac_behavior( subu )
{
  dbg_printf("subu r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] - RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction slt behavior method.
void ac_behavior( slt )
{	
  dbg_printf("slt r%d, r%d, r%d\n", rd, rs, rt);
  // Set the RD if RS< RT
  if( (ac_Sword) RB[rs] < (ac_Sword) RB[rt] )
    RB[rd] = 1;
  // Else reset RD
  else
    RB[rd] = 0;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction sltu behavior method.
void ac_behavior( sltu )
{
  dbg_printf("sltu r%d, r%d, r%d\n", rd, rs, rt);
  // Set the RD if RS < RT
  if( RB[rs] < RB[rt] )
    RB[rd] = 1;
  // Else reset RD
  else
    RB[rd] = 0;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction instr_and behavior method.
void ac_behavior( instr_and )
{
  dbg_printf("instr_and r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] & RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction instr_or behavior method.
void ac_behavior( instr_or )
{
  dbg_printf("instr_or r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] | RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction instr_xor behavior method.
void ac_behavior( instr_xor )
{
  dbg_printf("instr_xor r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] ^ RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction instr_nor behavior method.
void ac_behavior( instr_nor )
{
  dbg_printf("nor r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = ~(RB[rs] | RB[rt]);
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction nop behavior method.
void ac_behavior( nop )
{  
  dbg_printf("nop\n");
};

//!Instruction sll behavior method.
void ac_behavior( sll )
{  
  dbg_printf("sll r%d, r%d, %d\n", rd, rs, shamt);
  RB[rd] = RB[rt] << shamt;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction srl behavior method.
void ac_behavior( srl )
{
  dbg_printf("srl r%d, r%d, %d\n", rd, rs, shamt);
  RB[rd] = RB[rt] >> shamt;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction sra behavior method.
void ac_behavior( sra )
{
  dbg_printf("sra r%d, r%d, %d\n", rd, rs, shamt);
  RB[rd] = (ac_Sword) RB[rt] >> shamt;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction sllv behavior method.
void ac_behavior( sllv )
{
  dbg_printf("sllv r%d, r%d, r%d\n", rd, rt, rs);
  RB[rd] = RB[rt] << (RB[rs] & 0x1F);
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction srlv behavior method.
void ac_behavior( srlv )
{
  dbg_printf("srlv r%d, r%d, r%d\n", rd, rt, rs);
  RB[rd] = RB[rt] >> (RB[rs] & 0x1F);
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction srav behavior method.
void ac_behavior( srav )
{
  dbg_printf("srav r%d, r%d, r%d\n", rd, rt, rs);
  RB[rd] = (ac_Sword) RB[rt] >> (RB[rs] & 0x1F);
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction mult behavior method.
void ac_behavior( mult )
{
  dbg_printf("mult r%d, r%d\n", rs, rt);

  long long result;
  int half_result;

  result = (ac_Sword) RB[rs];
  result *= (ac_Sword) RB[rt];

  half_result = (result & 0xFFFFFFFF);
  // Register LO receives 32 less significant bits
  lo = half_result;

  half_result = ((result >> 32) & 0xFFFFFFFF);
  // Register HI receives 32 most significant bits
  hi = half_result ;

  dbg_printf("Result = %#llx\n", result);
};

//!Instruction multu behavior method.
void ac_behavior( multu )
{
  dbg_printf("multu r%d, r%d\n", rs, rt);

  unsigned long long result;
  unsigned int half_result;

  result  = RB[rs];
  result *= RB[rt];

  half_result = (result & 0xFFFFFFFF);
  // Register LO receives 32 less significant bits
  lo = half_result;

  half_result = ((result>>32) & 0xFFFFFFFF);
  // Register HI receives 32 most significant bits
  hi = half_result ;

  dbg_printf("Result = %#llx\n", result);
};

//!Instruction div behavior method.
void ac_behavior( div )
{
  dbg_printf("div r%d, r%d\n", rs, rt);
  // Register LO receives quotient
  lo = (ac_Sword) RB[rs] / (ac_Sword) RB[rt];
  // Register HI receives remainder
  hi = (ac_Sword) RB[rs] % (ac_Sword) RB[rt];
};

//!Instruction divu behavior method.
void ac_behavior( divu )
{
  dbg_printf("divu r%d, r%d\n", rs, rt);
  // Register LO receives quotient
  lo = RB[rs] / RB[rt];
  // Register HI receives remainder
  hi = RB[rs] % RB[rt];
};

//!Instruction mfhi behavior method.
void ac_behavior( mfhi )
{
  dbg_printf("mfhi r%d\n", rd);
  RB[rd] = hi;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction mthi behavior method.
void ac_behavior( mthi )
{
  dbg_printf("mthi r%d\n", rs);
  hi = RB[rs];
  dbg_printf("Result = %#x\n", (int)hi);
};

//!Instruction mflo behavior method.
void ac_behavior( mflo )
{
  dbg_printf("mflo r%d\n", rd);
  RB[rd] = lo;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction mtlo behavior method.
void ac_behavior( mtlo )
{
  dbg_printf("mtlo r%d\n", rs);
  lo = RB[rs];
  dbg_printf("Result = %#x\n", (int)lo);
};

//!Instruction j behavior method.
void ac_behavior( j )
{
  dbg_printf("j %d\n", addr);
  addr = addr << 2;
#ifndef NO_NEED_PC_UPDATE
  npc =  (ac_pc & 0xF0000000) | addr;
#endif 
  dbg_printf("Target = %#x\n", (ac_pc & 0xF0000000) | addr );
};

//!Instruction jal behavior method.
void ac_behavior( jal )
{
  dbg_printf("jal %d\n", addr);
  // Save the value of PC + 8 (return address) in $ra ($31) and
  // jump to the address given by PC(31...28)||(addr<<2)
  // It must also flush the instructions that were loaded into the pipeline
  RB[Ra] = ac_pc+4; //ac_pc is pc+4, we need pc+8
	
  addr = addr << 2;
#ifndef NO_NEED_PC_UPDATE
  npc = (ac_pc & 0xF0000000) | addr;
#endif 
	
  dbg_printf("Target = %#x\n", (ac_pc & 0xF0000000) | addr );
  dbg_printf("Return = %#x\n", ac_pc+4);
};

//!Instruction jr behavior method.
void ac_behavior( jr )
{
  dbg_printf("jr r%d\n", rs);
  // Jump to the address stored on the register reg[RS]
  // It must also flush the instructions that were loaded into the pipeline
#ifndef NO_NEED_PC_UPDATE
  npc = RB[rs], 1;
#endif 
  dbg_printf("Target = %#x\n", RB[rs]);
};

//!Instruction jalr behavior method.
void ac_behavior( jalr )
{
  dbg_printf("jalr r%d, r%d\n", rd, rs);
  // Save the value of PC + 8(return address) in rd and
  // jump to the address given by [rs]

#ifndef NO_NEED_PC_UPDATE
  npc = RB[rs], 1;
#endif 
  dbg_printf("Target = %#x\n", RB[rs]);

  if( rd == 0 )  //If rd is not defined use default
    rd = Ra;
  RB[rd] = ac_pc+4;
  dbg_printf("Return = %#x\n", ac_pc+4);
};

//!Instruction beq behavior method.
void ac_behavior( beq )
{
  dbg_printf("beq r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  if( RB[rs] == RB[rt] ){
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif
    if (predict == 0){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 1;
    }

    branchStalls += 3;
    branchStalls7 += 4;
    if (bp[(int)ac_pc] < 2){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      ++bp[(int)ac_pc];
    } else if (bp[(int)ac_pc] == 2){
      ++bp[(int)ac_pc];
    }
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  } else {
    if (predict == 1){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 0;
    }
    if (bp[(int)ac_pc] > 1){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]--;
    } else if (bp[(int)ac_pc] == 1){
      bp[(int)ac_pc]--;
    }
  }
};

//!Instruction bne behavior method.
void ac_behavior( bne )
{	
  dbg_printf("bne r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  if( RB[rs] != RB[rt] ){
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif
    if (predict == 0){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 1;
    }
    branchStalls += 3;
    branchStalls7 += 4;
     if (bp[(int)ac_pc] < 2){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]++;
    } else if (bp[(int)ac_pc] == 2){
      bp[(int)ac_pc]++;
    }
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  } else {
    if (predict == 1){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 0;
    }
    if (bp[(int)ac_pc] > 1){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]--;
    } else if (bp[(int)ac_pc] == 1){
      bp[(int)ac_pc]--;
    }
  }	
};

//!Instruction blez behavior method.
void ac_behavior( blez )
{
  dbg_printf("blez r%d, %d\n", rs, imm & 0xFFFF);
  if( (RB[rs] == 0 ) || (RB[rs]&0x80000000 ) ){
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2), 1;
#endif
    if (predict == 0){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 1;
    }
    branchStalls += 3;
    branchStalls7 += 4;
     if (bp[(int)ac_pc] < 2){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]++;
    } else if (bp[(int)ac_pc] == 2){
      bp[(int)ac_pc]++;
    }
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  } else {
    if (predict == 1){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 0;
    }
    if (bp[(int)ac_pc] > 1){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]--;
    } else if (bp[(int)ac_pc] == 1){
      bp[(int)ac_pc]--;
    }
  }	
};

//!Instruction bgtz behavior method.
void ac_behavior( bgtz )
{
  dbg_printf("bgtz r%d, %d\n", rs, imm & 0xFFFF);
  if( !(RB[rs] & 0x80000000) && (RB[rs]!=0) ){
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif
    if (predict == 0){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 1;
    }
    branchStalls += 3;
    branchStalls7 += 4;
     if (bp[(int)ac_pc] < 2){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]++;
    } else if (bp[(int)ac_pc] == 2){
      bp[(int)ac_pc]++;
    }
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  } else {
    if (predict == 1){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 0;
    }
    if (bp[(int)ac_pc] > 1){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]--;
    } else if (bp[(int)ac_pc] == 1){
      bp[(int)ac_pc]--;
    }
  }	
};

//!Instruction bltz behavior method.
void ac_behavior( bltz )
{
  dbg_printf("bltz r%d, %d\n", rs, imm & 0xFFFF);
  if( RB[rs] & 0x80000000 ){
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif
    if (predict == 0){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 1;
    }
    branchStalls += 3;
    branchStalls7 += 4;
     if (bp[(int)ac_pc] < 2){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]++;
    } else if (bp[(int)ac_pc] == 2){
      bp[(int)ac_pc]++;
    }
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  } else {
    if (predict == 1){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 0;
    }
    if (bp[(int)ac_pc] > 1){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]--;
    } else if (bp[(int)ac_pc] == 1){
      bp[(int)ac_pc]--;
    }
  }	
};

//!Instruction bgez behavior method.
void ac_behavior( bgez )
{
  dbg_printf("bgez r%d, %d\n", rs, imm & 0xFFFF);
  if( !(RB[rs] & 0x80000000) ){
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif
    if (predict == 0){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 1;
    }
    branchStalls += 3;
    branchStalls7 += 4;
     if (bp[(int)ac_pc] < 2){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]++;
    } else if (bp[(int)ac_pc] == 2){
      bp[(int)ac_pc]++;
    }
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  } else {
    if (predict == 1){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 0;
    }
    if (bp[(int)ac_pc] > 1){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]--;
    } else if (bp[(int)ac_pc] == 1){
      bp[(int)ac_pc]--;
    }
  }	
};

//!Instruction bltzal behavior method.
void ac_behavior( bltzal )
{
  dbg_printf("bltzal r%d, %d\n", rs, imm & 0xFFFF);
  RB[Ra] = ac_pc+4; //ac_pc is pc+4, we need pc+8
  if( RB[rs] & 0x80000000 ){
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif
    if (predict == 0){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 1;
    }
    branchStalls += 3;
    branchStalls7 += 4;
     if (bp[(int)ac_pc] < 2){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]++;
    } else if (bp[(int)ac_pc] == 2){
      bp[(int)ac_pc]++;
    }
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  } else {
    if (predict == 1){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 0;
    }
    if (bp[(int)ac_pc] > 1){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]--;
    } else if (bp[(int)ac_pc] == 1){
      bp[(int)ac_pc]--;
    }
  }	
  dbg_printf("Return = %#x\n", ac_pc+4);
};

//!Instruction bgezal behavior method.
void ac_behavior( bgezal )
{
  dbg_printf("bgezal r%d, %d\n", rs, imm & 0xFFFF);
  RB[Ra] = ac_pc+4; //ac_pc is pc+4, we need pc+8
  if( !(RB[rs] & 0x80000000) ){
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif
    if (predict == 0){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 1;
    }
    branchStalls += 3;
    branchStalls7 += 4;
     if (bp[(int)ac_pc] < 2){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]++;
    } else if (bp[(int)ac_pc] == 2){
      bp[(int)ac_pc]++;
    }
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  } else {
    if (predict == 1){
        dynamicStalls2 += 3;
        dynamicStalls27 += 4;
        bhist[(int)ac_pc][bp[(int)ac_pc]] = 0;
    }
    if (bp[(int)ac_pc] > 1){
      dynamicStalls += 3;
      dynamicStalls7 += 4;
      bp[(int)ac_pc]--;
    } else if (bp[(int)ac_pc] == 1){
      bp[(int)ac_pc]--;
    }
  }	
  dbg_printf("Return = %#x\n", ac_pc+4);
};

//!Instruction sys_call behavior method.
void ac_behavior( sys_call )
{
  dbg_printf("syscall\n");
  stop();
}

//!Instruction instr_break behavior method.
void ac_behavior( instr_break )
{
  fprintf(stderr, "instr_break behavior not implemented.\n"); 
  exit(EXIT_FAILURE);
}
