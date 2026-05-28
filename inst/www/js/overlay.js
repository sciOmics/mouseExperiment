// Mouse background overlay and diagnostics tab disabling
// Extracted from inline JS in app.R (issue 3.3)

$(document).ready(function() {
  var layer = $('<div id="mouse-bg-layer">' +
    '<img id="mouse-bg-img1" class="mouse-bg-img" src="images/mouse_outline1.png" alt="">' +
    '<img id="mouse-bg-img2" class="mouse-bg-img" src="images/mouse_outline2.png" alt="">' +
    '</div>');
  $('body').append(layer);

  // Disable Diagnostics tab — pending re-implementation
  function disableDiagnosticsTab() {
    var $diagLi = $('#tg_results_tabs li').filter(function() {
      return $(this).find('a[data-value="Diagnostics"]').length > 0;
    });
    if ($diagLi.length) {
      $diagLi.addClass('disabled').attr('title', 'Diagnostics — coming in a future release');
      $diagLi.find('a').css({'pointer-events': 'none', 'cursor': 'not-allowed'});
    }
  }
  $(document).one('shiny:sessioninitialized', disableDiagnosticsTab);
});

window.showAnalysisProgress = function(message, detail) {
  var el = document.getElementById('analysis-progress-overlay');
  if (!el) {
    el = document.createElement('div');
    el.id = 'analysis-progress-overlay';

    // Inline styles for all layout-critical properties so CSS cascade
    // cannot interfere with the overlay or card visibility.
    el.setAttribute('style',
      'position:fixed;top:0;left:0;width:100%;height:100%;' +
      'background:rgba(0,0,0,0.55);z-index:99990;' +
      'display:flex;align-items:center;justify-content:center;'
    );

    var isDark = document.body.classList.contains('dark-mode');
    var cardBg  = isDark ? '#2a2a2a' : '#ffffff';
    var msgCol  = isDark ? '#e0e0e0' : '#333333';
    var detCol  = isDark ? '#aaaaaa' : '#666666';
    var barBg   = isDark ? '#444444' : '#e9ecef';

    el.innerHTML =
      '<div id="analysis-progress-card" style="' +
          'background:' + cardBg + ';' +
          'border-radius:12px;' +
          'box-shadow:0 8px 32px rgba(0,0,0,0.45);' +
          'padding:36px 48px;' +
          'min-width:420px;max-width:560px;' +
          'text-align:center;">' +
        '<div id="analysis-progress-message" style="' +
            'font-size:1.15em;font-weight:600;' +
            'color:' + msgCol + ';margin-bottom:18px;"></div>' +
        '<div id="analysis-progress-bar-wrap" style="' +
            'height:6px;background:' + barBg + ';' +
            'border-radius:3px;overflow:hidden;' +
            'margin-bottom:14px;position:relative;">' +
          '<div id="analysis-progress-bar" style="' +
              'height:100%;width:5%;' +
              'background-color:var(--primary-color,#2d6a4f);' +
              'border-radius:3px;position:relative;overflow:hidden;' +
              'animation:progress-shimmer 1.6s ease-in-out infinite;"></div>' +
        '</div>' +
        '<div id="analysis-progress-detail" style="' +
            'font-size:0.9em;color:' + detCol + ';"></div>' +
      '</div>';

    document.body.appendChild(el);
  }

  document.getElementById('analysis-progress-message').textContent = message || '';
  document.getElementById('analysis-progress-detail').textContent  = detail  || '';
};

window.hideAnalysisProgress = function() {
  var el = document.getElementById('analysis-progress-overlay');
  if (el) el.style.display = 'none';
};
