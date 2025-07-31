#  run_BSA.sh
#  
#  Copyright 2021 WangPF <wangpf0608@126.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
#!/bin/bash
# Author:
#	WangPF
# Program:
# 	BSA 
# History:
# 	2021-02-20	First release  
# 	2023-03-31	Second release
## 加载配置文件
. ./.conf

## call variation
cd ${work_dir}/script
sh call_vari.sh
sh bgAnalysis.sh
sh statistics.sh
