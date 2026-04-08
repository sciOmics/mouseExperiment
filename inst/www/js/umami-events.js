/**
 * Umami Event Tracking Helper
 * Provides utility functions for tracking user interactions with Umami analytics
 */

// Wait for Umami to initialize before defining tracking functions
window.umamiReady = function() {
  return typeof window.umami !== 'undefined' && window.umami.track;
};

// Ensure we wait for Umami to be available (max 5 seconds)
window.umamiInit = new Promise(function(resolve) {
  var attempts = 0;
  var maxAttempts = 50; // 5 seconds with 100ms checks
  var checkUmami = setInterval(function() {
    attempts++;
    if (window.umamiReady() || attempts >= maxAttempts) {
      clearInterval(checkUmami);
      resolve(window.umamiReady());
    }
  }, 100);
});

// Tracking function wrapper (safe even if Umami is unavailable)
window.trackEvent = function(eventName, properties) {
  if (window.umamiReady()) {
    try {
      window.umami.track(eventName, properties || {});
      console.log('[Umami]', eventName, properties || {});
    } catch (e) {
      console.error('[Umami] Error tracking event:', eventName, e);
    }
  }
};

// Specialized tracking functions

/**
 * Track data file upload
 * @param {string} format - 'csv' or 'rda'
 * @param {number} rows - number of rows in dataset
 * @param {number} cols - number of columns
 */
window.trackDataUpload = function(format, rows, cols) {
  window.trackEvent('data_upload', {
    format: format,
    rows: rows,
    columns: cols,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track column mapping completion
 * @param {string} treatment - treatment column name
 * @param {string} time - time column name
 * @param {string} volume - volume column name
 * @param {boolean} survivalMapped - whether survival event column is mapped
 */
window.trackColumnMapping = function(treatment, time, volume, survivalMapped) {
  window.trackEvent('column_mapping_complete', {
    treatment_col: treatment,
    time_col: time,
    volume_col: volume,
    survival_event_mapped: survivalMapped,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track tab navigation
 * @param {string} tabName - name of tab being navigated to
 */
window.trackTabNavigation = function(tabName) {
  window.trackEvent('tab_navigate', {
    tab: tabName,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track analysis start
 * @param {string} analysisType - 'tumor_growth', 'survival', etc.
 * @param {string} method - 'lme4', 'auc', 'cox', etc.
 * @param {string} referenceGroup - reference group name
 */
window.trackAnalysisStart = function(analysisType, method, referenceGroup) {
  window.analysisStartTime = Date.now(); // Store start time for duration calculation
  window.trackEvent('analysis_run', {
    analysis_type: analysisType,
    method: method,
    reference_group: referenceGroup,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track analysis completion with runtime
 * @param {string} analysisType - 'tumor_growth', 'survival', etc.
 * @param {boolean} success - whether analysis completed successfully
 * @param {string} errorMessage - error message if failed (optional)
 * @param {number} durationSeconds - optional override for duration (otherwise uses start time)
 */
window.trackAnalysisComplete = function(analysisType, success, errorMessage, durationSeconds) {
  var duration = durationSeconds;
  if (!duration && window.analysisStartTime) {
    duration = (Date.now() - window.analysisStartTime) / 1000;
    window.analysisStartTime = null;
  }
  
  var properties = {
    analysis_type: analysisType,
    success: success,
    duration_seconds: duration || 0,
    timestamp: new Date().toISOString()
  };
  
  if (errorMessage) {
    properties.error = errorMessage;
  }
  
  window.trackEvent('analysis_complete', properties);
};

/**
 * Track plot generation
 * @param {string} analysisType - 'tumor_growth', 'survival', etc.
 * @param {string} plotType - 'growth_curve', 'km_curve', 'forest_plot', etc.
 */
window.trackPlotGenerated = function(analysisType, plotType) {
  window.trackEvent('plot_generated', {
    analysis_type: analysisType,
    plot_type: plotType,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track download action
 * @param {string} analysisType - 'tumor_growth', 'survival', etc.
 * @param {string} contentType - 'plot', 'table', 'results', etc.
 * @param {string} format - 'png', 'pdf', 'csv', 'txt', etc.
 */
window.trackDownload = function(analysisType, contentType, format) {
  window.trackEvent('download_results', {
    analysis_type: analysisType,
    content_type: contentType,
    format: format,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track error condition
 * @param {string} errorType - 'validation', 'analysis', 'file_read', etc.
 * @param {string} errorMessage - description of error
 * @param {string} context - module/section where error occurred
 */
window.trackError = function(errorType, errorMessage, context) {
  window.trackEvent('error_encountered', {
    error_type: errorType,
    error_message: errorMessage.substring(0, 200), // Truncate to 200 chars
    context: context,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track reference group selection
 * @param {string} analysisType - 'tumor_growth' or 'survival'
 * @param {string} referenceGroup - selected reference group value
 */
window.trackReferenceGroupSelected = function(analysisType, referenceGroup) {
  window.trackEvent('reference_group_selected', {
    analysis_type: analysisType,
    reference_group: referenceGroup,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track demo data loading
 * @param {string} demoName - name of demo dataset
 */
window.trackDemoDataLoaded = function(demoName) {
  window.trackEvent('demo_data_loaded', {
    demo_name: demoName,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track feature access
 * @param {string} featureName - name of feature being accessed
 * @param {string} status - 'enabled', 'disabled', 'coming_soon', etc.
 */
window.trackFeatureAccess = function(featureName, status) {
  window.trackEvent('feature_access', {
    feature: featureName,
    status: status,
    timestamp: new Date().toISOString()
  });
};

console.log('[Umami Events] Helper functions loaded. Waiting for Umami script to initialize...');
