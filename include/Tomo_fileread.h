#pragma once
#ifndef TOMO_FILEREAD_H_
#define TOMO_FILEREAD_H_

#include "CommonHeaders.h"
#include "Tomo_srs2D.h"
#include "Tomo_fsm2D.h"
#include "Tomo_inverse2D.h"
#include "Tomo_mesh2D.h"

void TomostartupRead(string filename, int &raytype, int& fwd_only, string& Fwd_Setting_File, string& Inverse_Setting_File,
	string& Mesh_Setting_File, string& Observed_data_File);

#endif