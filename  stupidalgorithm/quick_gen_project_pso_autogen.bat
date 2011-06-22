@echo off
set script_type=autogen
set cwd=%~pd0
set toolkit_path=D:\Program Files\exDev\exVim\toolkit
set lang_type=c cpp c# java shader python vim uc matlab wiki ini make sh batch debug qt swig
set vimfiles_path=.vimfiles.pso
set file_filter=*.C *.C++ *.CC *.CP *.CPP *.CXX *.H *.H++ *.HH *.HP *.HPP *.HXX *.INL *.IPP *.CS *.JAVA *.HLSL *.VSH *.PSH *.FX *.FXH *.CG *.SHD *.GLSL *.PY *.PYW *.PYX *.PXD *.VIM *.UC *.M *.WIKI *.INI *.CFG *.MAK *.MK *.MAKEFILE *.SH *.BSH *.BASH *.KSH *.ZSH *.BAT *.LOG *.ERR *.EXE *.QRC *.PRO *.PRI *.I *.SWG 
set file_filter_pattern="\\.c$|\\.C$|\\.c++$|\\.cc$|\\.cp$|\\.cpp$|\\.cxx$|\\.h$|\\.H$|\\.h++$|\\.hh$|\\.hp$|\\.hpp$|\\.hxx$|\\.inl$|\\.ipp$|\\.cs$|\\.java$|\\.hlsl$|\\.vsh$|\\.psh$|\\.fx$|\\.fxh$|\\.cg$|\\.shd$|\\.glsl$|\\.py$|\\.pyw$|\\.pyx$|\\.pxd$|\\.vim$|\\.uc$|\\.m$|\\.wiki$|\\.ini$|\\.cfg$|\\.mak$|\\.mk$|\\.Makefile$|\\.makefile$|\\.sh$|\\.SH$|\\.bsh$|\\.bash$|\\.ksh$|\\.zsh$|\\.bat$|\\.log$|\\.err$|\\.exe$|\\.qrc$|\\.pro$|\\.pri$|\\.i$|\\.swg$"
set cscope_file_filter=*.C *.C++ *.CC *.CP *.CPP *.CXX *.H *.H++ *.HH *.HP *.HPP *.HXX *.INL *.IPP *.HLSL *.VSH *.PSH *.FX *.FXH *.CG *.SHD *.GLSL 
set cscope_file_filter_pattern="\\.c$|\\.C$|\\.c++$|\\.cc$|\\.cp$|\\.cpp$|\\.cxx$|\\.h$|\\.H$|\\.h++$|\\.hh$|\\.hp$|\\.hpp$|\\.hxx$|\\.inl$|\\.ipp$|\\.hlsl$|\\.vsh$|\\.psh$|\\.fx$|\\.fxh$|\\.cg$|\\.shd$|\\.glsl$"
set dir_filter=
set support_filenamelist=true
set support_ctags=true
set support_symbol=true
set support_inherit=true
set support_cscope=true
set support_idutils=true
set ctags_cmd=ctags
set ctags_options= --c-kinds=+p --c++-kinds=+p --fields=+iaS --extra=+q --languages=c,c++,c#,java,python,vim,matlab,make,sh,c, --langmap=c:+.C,c++:+.H,c++:+.inl,c++:+.ipp,python:+.pyw,c:+.hlsl,c:+.vsh,c:+.psh,c:+.fx,c:+.fxh,c:+.cg,c:+.shd,c:+.glsl,
if exist .\%vimfiles_path%\quick_gen_project_pre_custom.bat (
    call .\%vimfiles_path%\quick_gen_project_pre_custom.bat
)
call "%toolkit_path%\quickgen\batch\quick_gen_project.bat" %1
if exist .\%vimfiles_path%\quick_gen_project_post_custom.bat (
    call .\%vimfiles_path%\quick_gen_project_post_custom.bat
)
echo on
