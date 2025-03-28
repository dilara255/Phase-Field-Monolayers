-- TODO: review usage of pchheader

--SOLUTION: PFM
workspace "PFM"

	startproject "controlCL"
	defines { "F_V2_API=", "F_AUX_API=", "F_CLIENTAPP", "PFM_API=" } -- using static libs


	toolset("clang")
	flags { "MultiProcessorCompile", "Verbose" }
	warnings "Extra"

	configurations {"Debug", "Release"}

	platforms { "x86_64" }

	filter "platforms:x86_64"
		architecture "x86_64"
		defines "SYS_ARCH=x86_64"
	filter {}

	filter "system:windows"
      defines "F_OS_WINDOWS"
	filter "system:linux"
      defines "F_OS_LINUX"
      linSharedLibDest_x86_64 = "/usr/local/lib64/"
	filter {}

	cfgDir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"
	guiCfgRelPath = "../res/guiCfg/imgui.ini"

	binDir = "bin/" .. cfgDir .. "/%{prj.name}"
	binIntDir = "bin-int/" .. cfgDir .. "/%{prj.name}"

	IncludeDir = {}
	IncludeDir["F_AUX"]   = "%{wks.location}/depend/fAux/include/"
	IncludeDir["F_VIZ2D"] = "%{wks.location}/depend/fViz2D/include/"
	IncludeDir["PFM_SIMUL"]   = "%{wks.location}/simulation/API"

	LibDir = {}
	LibDir["F_AUX"]   = ("%{wks.location}/depend/fAux/lib/" .. cfgDir)
	LibDir["F_VIZ2D"] = ("%{wks.location}/depend/fViz2D/lib/" .. cfgDir)
	LibDir["PFM_SIMUL"] = ("%{wks.location}/simulation/lib/" .. cfgDir)

--PROJECT: simulation
project "simulation"
	location "simulation"
	kind "Staticlib"
	language "C++"
	cppdialect "C++17"
	staticruntime "off"
	pic "on"

	undefines "F_CLIENTAPP"
	defines "F_PFM_SIMUL"
	defines "AS_BUILD_LIB"

	targetdir (binDir)
	objdir (binIntDir)

	files
	{
		"%{prj.name}/src/**.cpp",
		"%{prj.name}/src/**.h",
		"%{prj.name}/src/**.hpp",
		"%{prj.name}/include/**.hpp",
		"%{prj.name}/include/**.h",
		"%{prj.name}/API/**.h",
		"%{prj.name}/API/**.hpp"
	}

	includedirs {
		"%{prj.name}/include",
		"%{prj.name}/API",
		"%{IncludeDir.F_AUX}"
	}

	filter "architecture:x86_64"
		defines "X64"
	filter {}

	filter "system:windows"
		systemversion "latest"
		--buildoptions "/MT" --may cause override, should do inside filter
	filter {}

	filter "configurations:Debug"
		defines "AS_DEBUG"
		symbols "on"
	filter "configurations:Release"
		defines	"AS_RELEASE"
		optimize "on"
	filter {}

	postbuildcommands{ ("{COPYDIR} ../" .. binDir .. " %{LibDir.PFM_SIMUL}") }

--PROJECT: controlCL
project "controlCL"
	location "controlCL"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++17"
	staticruntime "off"

	defines "PFM_CONTROL_CL"
	defines "AS_BUILD_APP"

	links ("simulation")
	links ("%{LibDir.F_AUX}/fAux")

	targetdir (binDir)
	objdir (binIntDir)

	files
	{
		"%{prj.name}/src/**.cpp",
		"%{prj.name}/src/**.hpp",
		"%{prj.name}/src/**.h",
		"%{prj.name}/include/**.h",
		"%{prj.name}/include/**.hpp"
	}

	includedirs
	{
		"%{prj.name}/include",
		"%{IncludeDir.F_AUX}",
		"%{IncludeDir.PFM_SIMUL}"	
	}

	filter "architecture:x86_64"
		defines "X64"
	filter {}

	filter "system:windows"
		systemversion "latest"
		--buildoptions "/MT" --may cause override, should do inside filter
	filter {}

	filter "configurations:Debug"
		defines "AS_DEBUG"
		symbols "on"
	filter "configurations:Release"
		defines	"AS_RELEASE"
		optimize "on"
	filter {}

--PROJECT: controlGUI
project "controlGUI"
	location "controlGUI"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++17"
	staticruntime "off"

	defines "PFM_CONTROL_GUI"
	defines "AS_BUILD_APP"

	links ("simulation")
	links ("%{LibDir.F_AUX}/fAux")

	filter "system:linux"
		links { "glfw", "OpenGL" } -- TEST: IS THIS NEEDED HERE ON WIN?
	filter {}

	links ("%{LibDir.F_VIZ2D}/fViz2D")

	targetdir (binDir)
	objdir (binIntDir)

	files
	{
		"%{prj.name}/src/**.cpp",
		"%{prj.name}/src/**.hpp",
		"%{prj.name}/src/**.h",
		"%{prj.name}/include/**.h",
		"%{prj.name}/include/**.hpp"
	}

	includedirs
	{
		"%{prj.name}/include",
		"%{IncludeDir.F_AUX}",
		"%{IncludeDir.F_VIZ2D}",
		"%{IncludeDir.PFM_SIMUL}"	
	}

	filter "architecture:x86_64"
		defines "X64"
	filter {}

	filter "system:windows"
		systemversion "latest"
		--buildoptions "/MT" --may cause override, should do inside filter
	filter {}

	filter "configurations:Debug"
		defines "AS_DEBUG"
		symbols "on"
	filter "configurations:Release"
		defines	"AS_RELEASE"
		optimize "on"
	filter {}
	
	postbuildcommands{ ("{COPYFILE} " .. guiCfgRelPath .. " ../bin/" .. cfgDir .."/controlGUI/") }
