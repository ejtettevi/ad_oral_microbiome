#!/usr/bin/env python3
"""
Pathway-Tools Installation Helper Script
This script helps you install and configure Pathway-Tools for advanced microbiome analysis.
"""

import os
import sys
import subprocess
import platform

def check_system():
    """Check the current system and provide installation instructions"""
    print("Pathway-Tools Installation Helper")
    print("=" * 40)
    
    system = platform.system()
    print(f"Detected system: {system}")
    
    if system == "Darwin":  # macOS
        print("\nFor macOS:")
        print("1. Download Pathway-Tools from: https://bioinformatics.ai.sri.com/ptools/")
        print("2. Extract to /Applications/pathway-tools/ or /usr/local/pathway-tools/")
        print("3. Add to PATH by adding this line to ~/.zshrc:")
        print("   export PATH=\"/Applications/pathway-tools:$PATH\"")
        print("4. Reload shell: source ~/.zshrc")
        
    elif system == "Linux":
        print("\nFor Linux:")
        print("1. Download Pathway-Tools from: https://bioinformatics.ai.sri.com/ptools/")
        print("2. Extract to /usr/local/pathway-tools/ or ~/pathway-tools/")
        print("3. Add to PATH by adding this line to ~/.bashrc:")
        print("   export PATH=\"/usr/local/pathway-tools:$PATH\"")
        print("4. Reload shell: source ~/.bashrc")
        
    elif system == "Windows":
        print("\nFor Windows:")
        print("1. Download Pathway-Tools from: https://bioinformatics.ai.sri.com/ptools/")
        print("2. Install to C:\\pathway-tools\\")
        print("3. Add C:\\pathway-tools\\ to your system PATH")
        
    else:
        print(f"\nUnsupported system: {system}")
        print("Please check the Pathway-Tools documentation for your system.")

def check_installation():
    """Check if Pathway-Tools is properly installed"""
    print("\nChecking Pathway-Tools installation...")
    
    try:
        result = subprocess.run(['pathway-tools', '--version'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print("✓ Pathway-Tools found!")
            print(f"  Version: {result.stdout.strip()}")
            return True
        else:
            print("✗ Pathway-Tools command failed")
            return False
    except (subprocess.TimeoutExpired, FileNotFoundError):
        print("✗ Pathway-Tools command not found in PATH")
        return False

def check_python_api():
    """Check if Pathway-Tools Python API is available"""
    print("\nChecking Pathway-Tools Python API...")
    
    try:
        import ptools
        print("✓ Pathway-Tools Python API found!")
        return True
    except ImportError:
        print("✗ Pathway-Tools Python API not found")
        print("  Note: The Python API may not be available for all installations")
        print("  The script will still work with command-line interface")
        return False

def test_pathway_tools():
    """Test Pathway-Tools functionality"""
    print("\nTesting Pathway-Tools functionality...")
    
    if not check_installation():
        print("Cannot test - Pathway-Tools not installed")
        return False
    
    try:
        # Test basic command
        result = subprocess.run(['pathway-tools', '--help'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print("✓ Basic Pathway-Tools functionality working")
            return True
        else:
            print("✗ Basic Pathway-Tools functionality failed")
            return False
    except Exception as e:
        print(f"✗ Error testing Pathway-Tools: {e}")
        return False

def provide_installation_help():
    """Provide detailed installation help"""
    print("\n" + "=" * 50)
    print("DETAILED INSTALLATION INSTRUCTIONS")
    print("=" * 50)
    
    print("\n1. DOWNLOAD PATHWAY-TOOLS")
    print("   - Visit: https://bioinformatics.ai.sri.com/ptools/")
    print("   - Download the appropriate version for your system")
    print("   - You may need to register for a free account")
    
    print("\n2. INSTALL PATHWAY-TOOLS")
    print("   - Extract the downloaded file")
    print("   - Move to a permanent location (e.g., /Applications/ on macOS)")
    print("   - Make sure the pathway-tools executable is in the bin/ directory")
    
    print("\n3. ADD TO PATH")
    print("   - Find the directory containing the pathway-tools executable")
    print("   - Add this directory to your system PATH")
    print("   - Restart your terminal or reload your shell configuration")
    
    print("\n4. VERIFY INSTALLATION")
    print("   - Run: pathway-tools --version")
    print("   - You should see version information")
    
    print("\n5. (OPTIONAL) INSTALL PYTHON API")
    print("   - Some Pathway-Tools installations include a Python API")
    print("   - Try: pip install ptools")
    print("   - If not available, the script will use command-line interface")
    
    print("\n6. TEST WITH YOUR SCRIPT")
    print("   - Run your microbiome analysis script")
    print("   - It will automatically detect and use Pathway-Tools if available")

def main():
    """Main function"""
    print("Pathway-Tools Installation Helper")
    print("=" * 40)
    
    # Check system
    check_system()
    
    # Check current installation
    ptools_installed = check_installation()
    python_api_available = check_python_api()
    
    # Test functionality
    if ptools_installed:
        test_pathway_tools()
    
    # Provide help
    if not ptools_installed:
        print("\n" + "=" * 50)
        print("INSTALLATION REQUIRED")
        print("=" * 50)
        provide_installation_help()
    else:
        print("\n" + "=" * 50)
        print("INSTALLATION STATUS")
        print("=" * 50)
        print("✓ Pathway-Tools is installed and working")
        if python_api_available:
            print("✓ Python API is available")
            print("  Your script will use advanced Pathway-Tools features!")
        else:
            print("⚠ Python API not available")
            print("  Your script will use command-line interface")
            print("  This is still functional but with limited features")
    
    print("\n" + "=" * 50)
    print("NEXT STEPS")
    print("=" * 50)
    if ptools_installed:
        print("1. Run your microbiome analysis script")
        print("2. It will automatically use Pathway-Tools features")
        print("3. Check the output for advanced metabolic analysis results")
    else:
        print("1. Follow the installation instructions above")
        print("2. Install Pathway-Tools and add to PATH")
        print("3. Run this script again to verify installation")
        print("4. Then run your microbiome analysis script")

if __name__ == "__main__":
    main() 