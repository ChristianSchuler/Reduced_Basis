/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   constEq.h
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

//---------------------------------------------------------------------------
//..... LOCAL LEVEL INTEGRATION ALGORITHMS FOR CONSTITUTIVE EQUATIONS  ......
//---------------------------------------------------------------------------
#ifndef __constEq_h__
#define __constEq_h__

//---------------------------------------------------------------------------

struct Material_t;
struct Soft_t;
struct Controls;
struct SolVarDev;
struct SolVarBulk;
struct SolVarCell;
struct SolVarEdge;
struct PData;
struct JacRes;
struct Ph_trans_t;
struct DBMat;
struct Scaling;
//---------------------------------------------------------------------------

// constitutive equations evaluation context
struct ConstEqCtx
{
	// database parameters
	PetscInt     numPhases; 	// number phases
	Material_t  *phases;    	// phase parameters
	Soft_t      *soft;      	// material softening laws
	Ph_trans_t  *PhaseTrans;    // Phase transition laws
	DBMat       *dbm;
	Controls    *ctrl;      	// parameters and controls
	PData       *Pd;        	// phase diagram data
	Scaling     *scal;      	// scaling
	PetscScalar  dt;        	// time step
	PetscScalar  stats[3];  	// total number of [starts, successes, iterations]
	PetscScalar  avg_topo;  	// average surface topography

	// control volume parameters
	PetscScalar *phRat;  // phase ratios in the control volume
	SolVarDev   *svDev;  // deviatoric variables
	SolVarBulk  *svBulk; // volumetric variables
	PetscScalar  p;      // pressure
	PetscScalar  p_lith; // lithostatic pressure
	PetscScalar  p_pore; // pore pressure
	PetscScalar  T;      // temperature
	PetscScalar  DII;    // effective strain rate
	PetscScalar  Le;     // characteristic element size
	PetscScalar  depth;  // depth for depth-dependent density model

	// phase parameters
	PetscScalar  A_els;  // elasticity constant
	PetscScalar  A_dif;  // diffusion constant
	PetscScalar  A_max;  // upper bound constant
	PetscScalar  A_dis;  // dislocation constant
	PetscScalar  N_dis;  // dislocation exponent
	PetscScalar  A_prl;  // Peierls constant
	PetscScalar  N_prl;  // Peierls exponent
	PetscScalar  taupl;  // plastic yield stress

	// control volume results
	PetscScalar  eta;    // effective viscosity
	PetscScalar  eta_cr; // creep viscosity
	PetscScalar  DIIdif; // diffusion creep strain rate
	PetscScalar  DIIdis; // dislocation creep strain rate
	PetscScalar  DIIprl; // Peierls creep strain rate
	PetscScalar  DIIpl;  // plastic strain rate
	PetscScalar  yield;  // yield stress
};

//---------------------------------------------------------------------------
// setup evaluation context
PetscErrorCode setUpConstEq(ConstEqCtx *ctx, JacRes *jr);

// setup control volume parameters
PetscErrorCode setUpCtrlVol(
	ConstEqCtx  *ctx,    // context
	PetscScalar *phRat,  // phase ratios in the control volume
	SolVarDev   *svDev,  // deviatoric variables
	SolVarBulk  *svBulk, // volumetric variables
	PetscScalar  p,      // pressure
	PetscScalar  p_lith, // lithostatic pressure
	PetscScalar  p_pore, // pore pressure
	PetscScalar  T,      // temperature
	PetscScalar  DII,    // effective strain rate
	PetscScalar  z,      // z-coordinate of control volume
	PetscScalar  Le);    // characteristic element size

// setup phase parameters for deviatoric constitutive equation
PetscErrorCode setUpPhase(ConstEqCtx *ctx, PetscInt ID);

// evaluate deviatoric constitutive equations in control volume
PetscErrorCode devConstEq(ConstEqCtx *ctx);

// compute phase viscosities and strain rate partitioning
PetscErrorCode getPhaseVisc(ConstEqCtx *ctx, PetscInt ID);

// compute residual of the visco-elastic constitutive equation
PetscScalar getConsEqRes(PetscScalar eta, void *pctx);

// apply strain softening to a parameter (friction, cohesion)
PetscScalar applyStrainSoft(
		Soft_t      *soft, // material softening laws
		PetscInt     ID,   // softening law ID
		PetscScalar  APS,  // accumulated plastic strain
		PetscScalar  Le,   // characteristic element size
		PetscScalar  par); // softening parameter

// compute inverse elastic parameter in control volume
PetscScalar getI2Gdt(
		PetscInt     numPhases, // number phases
		Material_t  *phases,    // phase parameters
		PetscScalar *phRat,     // phase ratios in the control volume
		PetscScalar  dt);       // time step

// evaluate volumetric constitutive equations in control volume
PetscErrorCode volConstEq(ConstEqCtx *ctx);

// evaluate constitutive equations on the cell
PetscErrorCode cellConstEq(
		ConstEqCtx  *ctx,    // evaluation context
		SolVarCell  *svCell, // solution variables
		PetscScalar  dxx,    // effective normal strain rate components
		PetscScalar  dyy,    // ...
		PetscScalar  dzz,    // ...
		PetscScalar &sxx,    // Cauchy stress components
		PetscScalar &syy,    // ...
		PetscScalar &szz,    // ...
		PetscScalar &gres,   // volumetric residual
		PetscScalar &rho);   // effective density

// evaluate constitutive equations on the edge
PetscErrorCode edgeConstEq(
		ConstEqCtx  *ctx,    // evaluation context
		SolVarEdge  *svEdge, // solution variables
		PetscScalar  d,      // effective shear strain rate component
		PetscScalar &s);     // Cauchy stress component

// check convergence of constitutive equations
PetscErrorCode checkConvConstEq(ConstEqCtx *ctx);

//---------------------------------------------------------------------------
//.............................. PHASE DIAGRAM  .............................
//---------------------------------------------------------------------------

PetscErrorCode setDataPhaseDiagram(
		PData       *pd,
		PetscScalar  p,
		PetscScalar  T,
		char         pdn[]);

//---------------------------------------------------------------------------
#endif
