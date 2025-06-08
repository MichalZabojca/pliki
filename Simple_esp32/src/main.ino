/*
    Copyright (c) 2018 Gurpreet Bal https://github.com/ardyesp/ESPWebDAV
    Copyright (c) 2020 David Gauchard https://github.com/d-a-v/ESPWebDAV
    All rights reserved.

    Redistribution and use in source and binary forms, with or without modification,
    are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    3. The name of the author may not be used to endorse or promote products
      derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
    SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
    OF SUCH DAMAGE.

*/

/*  Using the WebDAV server

    From windows file explorer,
        \\ESPWebDAV.local\
        or http://ESPWebDAV.local/
        or Map Network Drive -> Connect to http://ESPWebDAV.local/
           subst w: \\ESPWebDAV.local/DavWWWRoot
    From macOS Finder > command-K > http://ESPWebDAV.local/
        (do not select anonymous to have write access)
    From macOS cmdline:
        mkdir -p /tmp/esp; mount_webdav -S -i -v esp32 http://ESPWebDAV.local/ /tmp/esp && echo OK
    From linux:
        use mount -t davs2 http://ESPWebDAV.local/ /mnt/
        use gio/gvfs/nautilus/YourFileExplorer http://ESPWebDAV.local/

    When running emulation on host (script ./run), use one of these URLs instead:
        http://ESPWebDAV.local:9080/
        http://127.0.0.1:9080/
        http://local-ip-address:9080/
        subst w: \\ESPWebDAV.local@9080/DavWWWRoot
*/
#include "esp_vfs_fat.h"
#include "sdmmc_cmd.h"
#include "driver/sdmmc_host.h"
#include <WiFi.h>
#include <ESPmDNS.h>
#include <SPI.h>
#include <SD.h>
#include <ESPWebDAV.h>
FS &gfs = SD;
#define HOSTNAME    "ESPWebDAV"
#ifndef STASSID
#define STASSID "PLAY_Swiatlowodowy_20CB"
#define STAPSK "8AfGdF8q5K@X"
#endif

#define FILESYSTEM SD
//#define FILESYSTEM SPIFFS

//WiFiServerSecure tcp(443);
WiFiServer tcp(80);

ESPWebDAV dav;

// ------------------------
void setup()
{
    SPI.begin(18, 32, 23, 5);
    delay(5000);
    Serial.begin(115200);
    Serial.print("MOSI: ");
    Serial.println(MOSI);
    Serial.print("MISO: ");
    Serial.println(MISO);
    Serial.print("SCK: ");
    Serial.println(SCK);
    Serial.print("SS: ");
    Serial.println(SS);  
    // ------------------------
    WiFi.persistent(false);
    WiFi.setHostname(HOSTNAME);
    WiFi.mode(WIFI_STA);
    Serial.begin(115200);
    WiFi.begin(STASSID, STAPSK);
    Serial.println("Connecting to " STASSID " ...");
    
    // Wait for connection
    while (WiFi.status() != WL_CONNECTED)
    {
        delay(500);
    }

    Serial.println("");
    Serial.print("Connected to "); Serial.println(STASSID);
    Serial.print("IP address: "); Serial.println(WiFi.localIP());
    Serial.print("RSSI: "); Serial.println(WiFi.RSSI());

    MDNS.begin(HOSTNAME);
    /*
    if (!FILESYSTEM.begin(5)) {
        Serial.println("Failed to mount SD");
        while (1) delay(1000);
    }
    tcp.begin();
    dav.begin(&tcp, &gfs);
    dav.setTransferStatusCallback([](const char* name, int percent, bool receive)
    {
        Serial.printf("%s: '%s': %d%%\n", receive ? "recv" : "send", name, percent);
    });

    Serial.println("WebDAV server started");
    */
    sdmmc_slot_config_t slot_config = SDMMC_SLOT_CONFIG_DEFAULT();
    sdmmc_card_t* card;
    const char mount_point[] = "/sdcard";
    slot_config.width = 1;
    esp_vfs_fat_sdmmc_mount_config_t mount_config = {
            .format_if_mount_failed = false,
            .max_files = 5,
            .allocation_unit_size = 16 * 1024
    };
    sdmmc_host_t host = SDMMC_HOST_DEFAULT();
    //host.flags = SDMMC_HOST_FLAG_4BIT;
    Serial.println("Mounting filesystem");
    esp_err_t ret = esp_vfs_fat_sdmmc_mount(mount_point, &host, &slot_config, &mount_config, &card);

    if (ret != ESP_OK) {
        if (ret == ESP_FAIL) {
            Serial.println("Failed to mount filesystem.");
        } 
        else {
            Serial.println("Failed to initialize the card (%s)."
                     "Make sure SD card lines have pull-up resistors in place.");
    /*
    #ifdef CONFIG_EXAMPLE_DEBUG_PIN_CONNECTIONS
            check_sd_card_pins(&config, pin_count);
    #endif
    */  
        }
        return;
    }
    Serial.println("Filesystem mounted");

    sdmmc_card_print_info(stdout, card);
    

}

void help()
{
    Serial.printf("interactive: F/ormat D/ir C/reateFile\n");

    Serial.printf("Heap stats: free heap: %u\n",
                  ESP.getHeapSize());
}

// ------------------------
void loop()
{
    dav.handleClient();

    int c = Serial.read();
    if (c > 0)
    {
        /*
        if (c == 'F')
        {
            Serial.println("formatting...");
            if (FILESYSTEM.format())
                Serial.println("Success");
            else
                Serial.println("Failure");
        }
        */
        if (c == 'D')
        {
            Serial.printf(">>>>>>>> dir /\n");
            dav.dir("/", &Serial);
            Serial.printf("<<<<<<<< dir\n");
        }
        else if (c == 'C')
        {
            auto f = gfs.open("/readme.md", "w");
            f.printf("hello\n");
            f.close();
        }
        else
            help();
    }

    // ------------------------
}
